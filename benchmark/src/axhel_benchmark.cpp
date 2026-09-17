// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/axhel.hpp"
#include "axhel/util/defines.hpp"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <linux/perf_event.h>
#include <sched.h>
#include <sys/prctl.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using unipi::axhel::EltwiseAddMod;
using unipi::axhel::EltwiseFMAMod;
using unipi::axhel::EltwiseMulMod;
using unipi::axhel::EltwiseReduceMod;
using unipi::axhel::EltwiseSubMod;
using unipi::axhel::NTT;

constexpr std::uint64_t kDefaultSeed = 0x415848454C42454EULL; // "AXHELBEN"

struct Config {
    std::vector<std::uint64_t> degrees{2048, 4096, 8192, 16384, 32768};
    std::vector<unsigned> bits{20, 25, 30, 40, 50, 52, 54, 56, 58, 60, 61, 62};
    std::vector<std::string> kernels;
    std::uint64_t repetitions{100};
    std::uint64_t warmup{10};
    // Number of consecutive public-kernel invocations measured within one PMU sample.
    // batch=1 preserves the original single-invocation measurement path.
    std::uint64_t batch{1};
    // PMU diagnostic profile: none | l1 | l2 | stalls.
    // none preserves the original single-cycle-counter benchmark path.
    std::string pmu_profile{"none"};
    // Optional diagnostic wall/reference timer. Disabled by default because
    // reading CNTVCT immediately around a very short kernel can perturb it.
    bool read_cntvct{false};
    std::uint64_t seed{kDefaultSeed};
    // Randomize only the order of (degree, qi_bits) scenarios. Operand generation
    // remains controlled by seed so datasets are unchanged across orderings.
    bool shuffle_scenarios{false};
    std::uint64_t shuffle_seed{0};
    int cpu{-2}; // -2 = auto first allowed, -1 = no pin, >=0 = requested CPU
    std::string output_dir{"benchmark-results"};
    std::string expect_backend;
    bool write_raw{true};
};

struct CounterReading {
    std::uint64_t raw_cycles{0};
    std::uint64_t time_enabled{0};
    std::uint64_t time_running{0};
    double cycles{0.0};
    bool scaled{false};

    // Diagnostic-only fields. They remain zero when --pmu-profile none.
    std::uint64_t cntvct_ticks{0};
    double instructions{0.0};
    double l1d_cache_rd{0.0};
    double l1d_cache_refill{0.0};
    double l2d_cache_rd{0.0};
    double l2d_cache_refill{0.0};
    double stall_backend{0.0};
    double stall_backend_mem{0.0};
    double stall_frontend{0.0};
};

class PerfCycles {
public:
    PerfCycles() {
        perf_event_attr attr{};
        attr.type = PERF_TYPE_HARDWARE;
        attr.size = sizeof(attr);
        attr.config = PERF_COUNT_HW_CPU_CYCLES;
        attr.disabled = 1;
        attr.exclude_kernel = 1;
        attr.exclude_hv = 1;
        attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;

        fd_ = static_cast<int>(::syscall(__NR_perf_event_open, &attr, 0, -1, -1, 0));
        if (fd_ < 0) {
            std::ostringstream oss;
            oss << "perf_event_open(PERF_COUNT_HW_CPU_CYCLES) failed: "
                << std::strerror(errno)
                << ". Check Linux perf permissions (kernel.perf_event_paranoid) and PMU access.";
            throw std::runtime_error(oss.str());
        }
    }

    PerfCycles(const PerfCycles&) = delete;
    PerfCycles& operator=(const PerfCycles&) = delete;

    ~PerfCycles() {
        if (fd_ >= 0) ::close(fd_);
    }

    void start() {
        if (::ioctl(fd_, PERF_EVENT_IOC_RESET, 0) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_RESET failed");
        }
        if (::ioctl(fd_, PERF_EVENT_IOC_ENABLE, 0) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_ENABLE failed");
        }
    }

    CounterReading stop() {
        if (::ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_DISABLE failed");
        }

        struct ReadData {
            std::uint64_t value;
            std::uint64_t time_enabled;
            std::uint64_t time_running;
        } data{};

        const auto n = ::read(fd_, &data, sizeof(data));
        if (n != static_cast<ssize_t>(sizeof(data))) {
            throw std::runtime_error("read(perf_event) failed");
        }
        // PERF_EVENT_IOC_RESET resets the event count, while time_enabled and
        // time_running can remain cumulative across enable/disable intervals.
        // Convert them to per-sample deltas before applying multiplex scaling.
        const std::uint64_t enabled_delta =
            data.time_enabled >= last_time_enabled_ ? data.time_enabled - last_time_enabled_ : data.time_enabled;
        const std::uint64_t running_delta =
            data.time_running >= last_time_running_ ? data.time_running - last_time_running_ : data.time_running;
        last_time_enabled_ = data.time_enabled;
        last_time_running_ = data.time_running;

        if (running_delta == 0) {
            throw std::runtime_error("hardware cycle counter was not scheduled during the sample");
        }

        CounterReading r;
        r.raw_cycles = data.value;
        r.time_enabled = enabled_delta;
        r.time_running = running_delta;
        r.scaled = running_delta != enabled_delta;
        if (r.scaled) {
            r.cycles = static_cast<double>(data.value) *
                       static_cast<double>(enabled_delta) /
                       static_cast<double>(running_delta);
        } else {
            r.cycles = static_cast<double>(data.value);
        }
        return r;
    }

private:
    int fd_{-1};
    std::uint64_t last_time_enabled_{0};
    std::uint64_t last_time_running_{0};
};

#if defined(__aarch64__)
static inline std::uint64_t read_cntvct_el0() {
    std::uint64_t value = 0;
    asm volatile("isb\n\tmrs %0, cntvct_el0" : "=r"(value) : : "memory");
    return value;
}

static inline std::uint64_t read_cntfrq_el0() {
    std::uint64_t value = 0;
    asm volatile("mrs %0, cntfrq_el0" : "=r"(value));
    return value;
}
#else
static inline std::uint64_t read_cntvct_el0() { return 0; }
static inline std::uint64_t read_cntfrq_el0() { return 0; }
#endif

unsigned read_arm_pmu_type() {
    const char* path = "/sys/bus/event_source/devices/armv8_pmuv3_0/type";
    std::ifstream in(path);
    unsigned type = 0;
    if (!(in >> type)) {
        throw std::runtime_error(
            std::string("cannot read Arm core PMU type from ") + path +
            ". This diagnostic mode expects the armv8_pmuv3_0 PMU exposed by perf.");
    }
    return type;
}

struct PmuEventSpec {
    const char* name;
    std::uint64_t config;
};

class PerfGroup {
public:
    explicit PerfGroup(const std::string& profile) : profile_(profile), pmu_type_(read_arm_pmu_type()) {
        // Architectural Arm PMUv3 event numbers plus the Neoverse-V2 event
        // reported by `perf list --details` on the target GH200 system.
        events_.push_back({"cycles", 0x11});
        events_.push_back({"instructions", 0x08});

        if (profile_ == "l1") {
            events_.push_back({"l1d_cache_rd", 0x40});
            events_.push_back({"l1d_cache_refill", 0x03});
        } else if (profile_ == "l2") {
            events_.push_back({"l2d_cache_rd", 0x50});
            events_.push_back({"l2d_cache_refill", 0x17});
        } else if (profile_ == "stalls") {
            events_.push_back({"stall_backend", 0x24});
            events_.push_back({"stall_backend_mem", 0x4005});
            events_.push_back({"stall_frontend", 0x23});
        } else {
            throw std::invalid_argument("unknown PMU diagnostic profile: " + profile_);
        }

        open_group();
    }

    PerfGroup(const PerfGroup&) = delete;
    PerfGroup& operator=(const PerfGroup&) = delete;

    ~PerfGroup() {
        for (auto it = fds_.rbegin(); it != fds_.rend(); ++it) {
            if (*it >= 0) ::close(*it);
        }
    }

    void start() {
        if (::ioctl(leader_fd_, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_RESET(group) failed");
        }
        if (::ioctl(leader_fd_, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_ENABLE(group) failed");
        }
    }

    CounterReading stop() {
        if (::ioctl(leader_fd_, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP) == -1) {
            throw std::runtime_error("PERF_EVENT_IOC_DISABLE(group) failed");
        }

        // PERF_FORMAT_GROUP | TOTAL_TIME_ENABLED | TOTAL_TIME_RUNNING | ID:
        // nr, time_enabled, time_running, then (value,id) for each event.
        std::vector<std::uint64_t> data(3 + 2 * events_.size(), 0);
        const std::size_t bytes = data.size() * sizeof(std::uint64_t);
        const auto got = ::read(leader_fd_, data.data(), bytes);
        if (got != static_cast<ssize_t>(bytes)) {
            throw std::runtime_error("read(perf PMU group) failed");
        }
        if (data[0] != events_.size()) {
            throw std::runtime_error("unexpected number of events returned by perf PMU group");
        }

        const std::uint64_t enabled_total = data[1];
        const std::uint64_t running_total = data[2];
        const std::uint64_t enabled_delta =
            enabled_total >= last_time_enabled_ ? enabled_total - last_time_enabled_ : enabled_total;
        const std::uint64_t running_delta =
            running_total >= last_time_running_ ? running_total - last_time_running_ : running_total;
        last_time_enabled_ = enabled_total;
        last_time_running_ = running_total;

        if (running_delta == 0) {
            throw std::runtime_error(
                "PMU diagnostic group was not scheduled; reduce the number of grouped events");
        }

        std::map<std::uint64_t, std::uint64_t> values_by_id;
        for (std::size_t i = 0; i < events_.size(); ++i) {
            const std::uint64_t value = data[3 + 2 * i];
            const std::uint64_t id = data[3 + 2 * i + 1];
            values_by_id[id] = value;
        }

        const bool scaled = running_delta != enabled_delta;
        const double scale = scaled
            ? static_cast<double>(enabled_delta) / static_cast<double>(running_delta)
            : 1.0;

        CounterReading r;
        r.time_enabled = enabled_delta;
        r.time_running = running_delta;
        r.scaled = scaled;

        for (std::size_t i = 0; i < events_.size(); ++i) {
            const auto pos = values_by_id.find(ids_[i]);
            if (pos == values_by_id.end()) {
                throw std::runtime_error("perf PMU group read returned an unknown/missing event id");
            }
            const std::uint64_t raw = pos->second;
            const double value = static_cast<double>(raw) * scale;
            const std::string name = events_[i].name;

            if (name == "cycles") {
                r.raw_cycles = raw;
                r.cycles = value;
            } else if (name == "instructions") {
                r.instructions = value;
            } else if (name == "l1d_cache_rd") {
                r.l1d_cache_rd = value;
            } else if (name == "l1d_cache_refill") {
                r.l1d_cache_refill = value;
            } else if (name == "l2d_cache_rd") {
                r.l2d_cache_rd = value;
            } else if (name == "l2d_cache_refill") {
                r.l2d_cache_refill = value;
            } else if (name == "stall_backend") {
                r.stall_backend = value;
            } else if (name == "stall_backend_mem") {
                r.stall_backend_mem = value;
            } else if (name == "stall_frontend") {
                r.stall_frontend = value;
            }
        }
        return r;
    }

private:
    void open_group() {
        int group_fd = -1;
        for (std::size_t i = 0; i < events_.size(); ++i) {
            perf_event_attr attr{};
            attr.type = pmu_type_;
            attr.size = sizeof(attr);
            attr.config = events_[i].config;
            attr.disabled = (i == 0) ? 1 : 0;
            attr.exclude_kernel = 1;
            attr.exclude_hv = 1;
            attr.pinned = (i == 0) ? 1 : 0;
            attr.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_TOTAL_TIME_ENABLED |
                               PERF_FORMAT_TOTAL_TIME_RUNNING | PERF_FORMAT_ID;

            const int fd = static_cast<int>(
                ::syscall(__NR_perf_event_open, &attr, 0, -1, group_fd, 0));
            if (fd < 0) {
                std::ostringstream oss;
                oss << "perf_event_open failed for PMU event " << events_[i].name
                    << " (config=0x" << std::hex << events_[i].config << std::dec << "): "
                    << std::strerror(errno);
                throw std::runtime_error(oss.str());
            }
            if (i == 0) {
                leader_fd_ = fd;
                group_fd = fd;
            }
            fds_.push_back(fd);

            std::uint64_t id = 0;
            if (::ioctl(fd, PERF_EVENT_IOC_ID, &id) == -1) {
                throw std::runtime_error("PERF_EVENT_IOC_ID failed for PMU group member");
            }
            ids_.push_back(id);
        }
    }

    std::string profile_;
    unsigned pmu_type_{0};
    int leader_fd_{-1};
    std::vector<PmuEventSpec> events_;
    std::vector<int> fds_;
    std::vector<std::uint64_t> ids_;
    std::uint64_t last_time_enabled_{0};
    std::uint64_t last_time_running_{0};
};

struct RawSample {
    std::uint64_t repetition{0};
    CounterReading counter;
};

struct Statistics {
    std::size_t samples_total{0};
    std::size_t samples_kept{0};
    std::size_t outliers_removed{0};
    double raw_min{0};
    double raw_max{0};
    double q1{0};
    double q3{0};
    double iqr{0};
    double lower_fence{0};
    double upper_fence{0};
    double min{0};
    double max{0};
    double mean{0};
    double median{0};
    double stddev{0};
    double cv_percent{0};
    bool any_scaled{false};
};

struct KernelMeta {
    std::string name;
    std::uint64_t mod_factor{0};
    std::uint64_t input_mod_factor{0};
    std::uint64_t output_mod_factor{0};
};

std::string backend_name() {
#ifdef AXHEL_BENCHMARK_STRICT_SCALAR
    return "ScalarStrict";
#elif defined(AXHEL_HAS_SVE)
    return "SVE";
#else
    return "Native";
#endif
}

std::string compiler_name() {
#if defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#else
    return "unknown";
#endif
}

std::string trim(std::string s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::vector<std::uint64_t> parse_u64_list(const std::string& s) {
    std::vector<std::uint64_t> out;
    for (const auto& item : split(s, ',')) out.push_back(std::stoull(item));
    return out;
}

std::vector<unsigned> parse_unsigned_list(const std::string& s) {
    std::vector<unsigned> out;
    for (const auto& item : split(s, ',')) out.push_back(static_cast<unsigned>(std::stoul(item)));
    return out;
}

const std::vector<std::string>& all_kernel_names() {
    static const std::vector<std::string> names{
        "add_vv", "add_vs",
        "sub_vv", "sub_vs",
        "mul_f1", "mul_f2", "mul_f4",
        "fma_f1", "fma_f2", "fma_f4", "fma_f8",
        "reduce_2to1", "reduce_4to1", "reduce_4to2", "reduce_barrett",
        "ntt_forward", "ntt_forward_lazy",
        "ntt_inverse", "ntt_inverse_lazy"
    };
    return names;
}

void print_help(const char* argv0) {
    std::cout
        << "AXHEL public-kernel HW-cycle benchmark\n\n"
        << "Usage: " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --repetitions N       measured repetitions per scenario (default 100)\n"
        << "  --warmup N            warm-up repetitions per scenario (default 10)\n"
        << "  --batch N             kernel invocations per PMU sample (default 1)\n"
        << "  --pmu-profile P       none|l1|l2|stalls (default none; diagnostic)\n"
        << "  --cntvct              also read CNTVCT_EL0 around each diagnostic invocation\n"
        << "  --degrees LIST        comma-separated polynomial degrees\n"
        << "  --bits LIST           comma-separated modulus bit widths\n"
        << "  --kernels LIST        comma-separated kernel names, or 'all'\n"
        << "  --cpu auto|off|ID     pin to first allowed CPU, no pin, or CPU ID\n"
        << "  --seed N              deterministic base seed for operand generation\n"
        << "  --shuffle-scenarios   randomize (degree, qi_bits) scenario order\n"
        << "  --shuffle-seed N      deterministic scenario-order seed\n"
        << "  --output-dir DIR      output directory\n"
        << "  --expect-backend B    require ScalarStrict, Native, or SVE (runner safety check)\n"
        << "  --no-raw              do not write raw_samples.csv\n"
        << "  --list-kernels        print supported kernel names\n"
        << "  --help                show this help\n";
}

Config parse_args(int argc, char** argv) {
    Config cfg;
    cfg.kernels = all_kernel_names();

    auto require_value = [&](int& i) -> std::string {
        if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value after ") + argv[i]);
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            print_help(argv[0]);
            std::exit(0);
        } else if (arg == "--list-kernels") {
            for (const auto& k : all_kernel_names()) std::cout << k << '\n';
            std::exit(0);
        } else if (arg == "--repetitions") {
            cfg.repetitions = std::stoull(require_value(i));
        } else if (arg == "--warmup") {
            cfg.warmup = std::stoull(require_value(i));
        } else if (arg == "--batch") {
            cfg.batch = std::stoull(require_value(i));
        } else if (arg == "--pmu-profile") {
            cfg.pmu_profile = require_value(i);
            if (cfg.pmu_profile != "none" && cfg.pmu_profile != "l1" &&
                cfg.pmu_profile != "l2" && cfg.pmu_profile != "stalls") {
                throw std::invalid_argument("--pmu-profile must be none, l1, l2, or stalls");
            }
        } else if (arg == "--cntvct") {
            cfg.read_cntvct = true;
        } else if (arg == "--degrees") {
            cfg.degrees = parse_u64_list(require_value(i));
        } else if (arg == "--bits") {
            cfg.bits = parse_unsigned_list(require_value(i));
        } else if (arg == "--kernels") {
            const auto value = require_value(i);
            cfg.kernels = (value == "all") ? all_kernel_names() : split(value, ',');
        } else if (arg == "--cpu") {
            const auto value = require_value(i);
            if (value == "auto") cfg.cpu = -2;
            else if (value == "off") cfg.cpu = -1;
            else cfg.cpu = std::stoi(value);
        } else if (arg == "--seed") {
            cfg.seed = std::stoull(require_value(i), nullptr, 0);
        } else if (arg == "--shuffle-scenarios") {
            cfg.shuffle_scenarios = true;
        } else if (arg == "--shuffle-seed") {
            cfg.shuffle_seed = std::stoull(require_value(i), nullptr, 0);
        } else if (arg == "--output-dir") {
            cfg.output_dir = require_value(i);
        } else if (arg == "--expect-backend") {
            cfg.expect_backend = require_value(i);
        } else if (arg == "--no-raw") {
            cfg.write_raw = false;
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    if (cfg.repetitions < 4) throw std::invalid_argument("--repetitions must be >= 4 for IQR statistics");
    if (cfg.batch == 0) throw std::invalid_argument("--batch must be >= 1");
    if (cfg.pmu_profile != "none" && cfg.batch != 1) {
        throw std::invalid_argument(
            "PMU diagnostic profiles require --batch 1 so every PMU row corresponds to one kernel invocation");
    }
    if (cfg.read_cntvct && cfg.pmu_profile == "none") {
        throw std::invalid_argument("--cntvct requires --pmu-profile l1, l2, or stalls");
    }
    if (cfg.degrees.empty() || cfg.bits.empty() || cfg.kernels.empty()) throw std::invalid_argument("empty benchmark matrix");

    for (const auto degree : cfg.degrees) {
        if (degree < 2 || (degree & (degree - 1)) != 0) {
            throw std::invalid_argument("all polynomial degrees must be powers of two >= 2");
        }
    }
    for (const auto bits : cfg.bits) {
        if (bits < 3 || bits > 62) throw std::invalid_argument("modulus bit widths must be in [3,62]");
    }
    for (const auto& k : cfg.kernels) {
        if (std::find(all_kernel_names().begin(), all_kernel_names().end(), k) == all_kernel_names().end()) {
            throw std::invalid_argument("unsupported kernel: " + k);
        }
    }

    // The NTT public kernels are in-place. Repeating them inside one measured
    // region would feed each invocation with the output of the previous one,
    // so batch>1 is intentionally restricted to element-wise kernels.
    if (cfg.batch > 1) {
        for (const auto& k : cfg.kernels) {
            if (k == "ntt_forward" || k == "ntt_forward_lazy" ||
                k == "ntt_inverse" || k == "ntt_inverse_lazy") {
                throw std::invalid_argument(
                    "--batch > 1 is supported only for element-wise kernels; "
                    "NTT kernels are in-place");
            }
        }
    }
    return cfg;
}

int pin_cpu(int request) {
    if (request == -1) return -1;

    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (::sched_getaffinity(0, sizeof(allowed), &allowed) != 0) {
        throw std::runtime_error(std::string("sched_getaffinity failed: ") + std::strerror(errno));
    }

    int cpu = request;
    if (request == -2) {
        cpu = -1;
        for (int i = 0; i < CPU_SETSIZE; ++i) {
            if (CPU_ISSET(i, &allowed)) {
                cpu = i;
                break;
            }
        }
        if (cpu < 0) throw std::runtime_error("no CPU available in process affinity mask");
    } else if (request < 0 || request >= CPU_SETSIZE || !CPU_ISSET(request, &allowed)) {
        throw std::runtime_error("requested CPU is not in the process affinity mask");
    }

    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (::sched_setaffinity(0, sizeof(set), &set) != 0) {
        throw std::runtime_error(std::string("sched_setaffinity failed: ") + std::strerror(errno));
    }
    return cpu;
}

std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b, std::uint64_t mod) {
    return static_cast<std::uint64_t>((static_cast<unsigned __int128>(a) * b) % mod);
}

std::uint64_t pow_mod(std::uint64_t a, std::uint64_t e, std::uint64_t mod) {
    std::uint64_t r = 1;
    while (e) {
        if (e & 1U) r = mul_mod(r, a, mod);
        a = mul_mod(a, a, mod);
        e >>= 1U;
    }
    return r;
}

bool is_prime64(std::uint64_t n) {
    if (n < 2) return false;
    for (const std::uint64_t p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    std::uint64_t d = n - 1;
    unsigned s = 0;
    while ((d & 1U) == 0) {
        d >>= 1U;
        ++s;
    }

    for (const std::uint64_t a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        std::uint64_t x = pow_mod(a % n, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (unsigned r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }
    return true;
}

std::uint64_t find_ntt_prime(std::uint64_t degree, unsigned bits) {
    const std::uint64_t step = 2 * degree;
    const std::uint64_t hi = (std::uint64_t{1} << bits) - 1;
    const std::uint64_t lo = std::uint64_t{1} << (bits - 1);
    std::uint64_t q = hi - ((hi - 1) % step);
    while (q >= lo && q > step) {
        if (is_prime64(q)) return q;
        if (q < lo + step) break;
        q -= step;
    }
    throw std::runtime_error("unable to find NTT-friendly prime for degree=" +
                             std::to_string(degree) + " bits=" + std::to_string(bits));
}

std::uint64_t find_primitive_2n_root(std::uint64_t degree, std::uint64_t q) {
    const std::uint64_t order = 2 * degree;
    const std::uint64_t exponent = (q - 1) / order;
    for (std::uint64_t a = 2; a < 100000; ++a) {
        const auto psi = pow_mod(a, exponent, q);
        if (psi != 1 && pow_mod(psi, degree, q) == q - 1 && pow_mod(psi, order, q) == 1) {
            return psi;
        }
    }
    throw std::runtime_error("unable to find primitive 2N-th root");
}

std::uint64_t bounded_random(std::mt19937_64& rng, std::uint64_t bound) {
    if (bound == 0) throw std::invalid_argument("bounded_random bound must be non-zero");
    const std::uint64_t threshold = static_cast<std::uint64_t>(-bound) % bound;
    for (;;) {
        const std::uint64_t x = rng();
        if (x >= threshold) return x % bound;
    }
}

std::uint64_t random_factor_value(std::mt19937_64& rng, std::uint64_t q, std::uint64_t factor) {
    const unsigned __int128 bound = static_cast<unsigned __int128>(q) * factor;
    const unsigned __int128 domain = static_cast<unsigned __int128>(std::numeric_limits<std::uint64_t>::max()) + 1;
    if (bound >= domain) return rng();
    return bounded_random(rng, static_cast<std::uint64_t>(bound));
}

std::vector<std::uint64_t> random_vector(std::size_t n, std::uint64_t q, std::uint64_t factor,
                                         std::mt19937_64& rng) {
    std::vector<std::uint64_t> v(n);
    for (auto& x : v) x = random_factor_value(rng, q, factor);
    return v;
}

std::uint64_t scenario_seed(std::uint64_t seed, std::uint64_t degree, unsigned bits, const std::string& kernel) {
    std::uint64_t x = seed ^ (degree + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U));
    x ^= static_cast<std::uint64_t>(bits) * 0xbf58476d1ce4e5b9ULL;
    for (const unsigned char c : kernel) {
        x ^= static_cast<std::uint64_t>(c) + 0x9e3779b97f4a7c15ULL + (x << 6U) + (x >> 2U);
    }
    return x;
}

double quantile_sorted(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) throw std::invalid_argument("quantile on empty data");
    if (sorted.size() == 1) return sorted.front();
    const double pos = p * static_cast<double>(sorted.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    const double frac = pos - static_cast<double>(lo);
    return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
}

Statistics compute_statistics(const std::vector<RawSample>& raw) {
    if (raw.empty()) throw std::invalid_argument("no samples");
    std::vector<double> values;
    values.reserve(raw.size());
    bool any_scaled = false;
    for (const auto& s : raw) {
        values.push_back(s.counter.cycles);
        any_scaled = any_scaled || s.counter.scaled;
    }
    std::sort(values.begin(), values.end());

    Statistics st;
    st.samples_total = values.size();
    st.raw_min = values.front();
    st.raw_max = values.back();
    st.q1 = quantile_sorted(values, 0.25);
    st.q3 = quantile_sorted(values, 0.75);
    st.iqr = st.q3 - st.q1;
    st.lower_fence = st.q1 - 1.5 * st.iqr;
    st.upper_fence = st.q3 + 1.5 * st.iqr;
    st.any_scaled = any_scaled;

    std::vector<double> kept;
    kept.reserve(values.size());
    std::copy_if(values.begin(), values.end(), std::back_inserter(kept), [&](double x) {
        return x >= st.lower_fence && x <= st.upper_fence;
    });
    if (kept.empty()) throw std::runtime_error("IQR filter removed every sample");

    st.samples_kept = kept.size();
    st.outliers_removed = st.samples_total - st.samples_kept;
    st.min = kept.front();
    st.max = kept.back();
    st.mean = std::accumulate(kept.begin(), kept.end(), 0.0) / static_cast<double>(kept.size());
    st.median = quantile_sorted(kept, 0.5);

    double sum_sq = 0.0;
    for (const double x : kept) {
        const double d = x - st.mean;
        sum_sq += d * d;
    }
    st.stddev = kept.size() > 1 ? std::sqrt(sum_sq / static_cast<double>(kept.size() - 1)) : 0.0;
    st.cv_percent = st.mean != 0.0 ? (100.0 * st.stddev / st.mean) : 0.0;
    return st;
}

template <typename Prepare, typename Invoke, typename Consume>
std::vector<RawSample> measure_kernel(std::uint64_t warmup, std::uint64_t repetitions,
                                      std::uint64_t batch, const std::string& pmu_profile,
                                      bool use_cntvct,
                                      Prepare&& prepare, Invoke&& invoke, Consume&& consume) {
    // Warm-up is outside every PMU measurement.
    for (std::uint64_t i = 0; i < warmup; ++i) {
        prepare();
        if (batch == 1) {
            invoke();
        } else {
            for (std::uint64_t j = 0; j < batch; ++j) invoke();
        }
        consume();
    }

    std::vector<RawSample> samples;
    samples.reserve(static_cast<std::size_t>(repetitions));

    if (pmu_profile == "none") {
        // IMPORTANT: this is the original measurement path. In particular,
        // --batch 1 executes no extra loop or CNTVCT read inside the PMU window.
        PerfCycles counter;
        for (std::uint64_t i = 0; i < repetitions; ++i) {
            prepare();

            CounterReading reading;
            if (batch == 1) {
                counter.start();
                invoke();
                reading = counter.stop();
            } else {
                counter.start();
                for (std::uint64_t j = 0; j < batch; ++j) invoke();
                reading = counter.stop();
                reading.cycles /= static_cast<double>(batch);
            }

            consume();
            samples.push_back({i, reading});
        }
        return samples;
    }

#if !defined(__aarch64__)
    (void)use_cntvct;
    throw std::runtime_error("--pmu-profile diagnostics require AArch64");
#else
    // Diagnostic path. batch==1 is enforced during option validation.
    PerfGroup counter(pmu_profile);
    for (std::uint64_t i = 0; i < repetitions; ++i) {
        prepare();

        CounterReading reading;
        if (use_cntvct) {
            counter.start();
            const std::uint64_t t0 = read_cntvct_el0();
            invoke();
            const std::uint64_t t1 = read_cntvct_el0();
            reading = counter.stop();
            reading.cntvct_ticks = t1 - t0;
        } else {
            // Preferred first diagnostic pass: no extra instruction is placed
            // between PMU enable and the public-kernel invocation.
            counter.start();
            invoke();
            reading = counter.stop();
        }

        consume();
        samples.push_back({i, reading});
    }
    return samples;
#endif
}

volatile std::uint64_t g_sink = 0;

int sve_vector_bits() {
#if defined(__aarch64__) && defined(PR_SVE_GET_VL) && defined(PR_SVE_VL_LEN_MASK)
    if (AXHEL_BENCHMARK_SVE_DETECTED) {
        const int value = ::prctl(PR_SVE_GET_VL);
        if (value >= 0) return (value & PR_SVE_VL_LEN_MASK) * 8;
    }
#endif
    return 0;
}

std::string kernel_release() {
    utsname u{};
    return ::uname(&u) == 0 ? std::string(u.release) : std::string("unknown");
}

std::string hostname() {
    utsname u{};
    return ::uname(&u) == 0 ? std::string(u.nodename) : std::string("unknown");
}

std::string read_cpu_model() {
    std::ifstream in("/proc/cpuinfo");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("model name", 0) == 0 || line.rfind("Processor", 0) == 0 || line.rfind("Hardware", 0) == 0) {
            const auto pos = line.find(':');
            if (pos != std::string::npos) return trim(line.substr(pos + 1));
        }
    }
    return "unknown";
}

void ensure_dir(const std::string& path) {
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    if (ec) throw std::runtime_error("unable to create output directory: " + path + ": " + ec.message());
}

void write_metadata(const Config& cfg, int pinned_cpu, const std::string& output_dir) {
    std::ofstream out(output_dir + "/metadata.txt");
    if (!out) throw std::runtime_error("cannot create metadata.txt");
    out << "backend=" << backend_name() << '\n';
#ifdef AXHEL_BENCHMARK_STRICT_SCALAR
    out << "strict_scalar=1\n";
#if defined(__clang__)
    out << "vectorization_policy=-fno-vectorize -fno-slp-vectorize\n";
#elif defined(__GNUC__)
    out << "vectorization_policy=-fno-tree-vectorize -fno-tree-loop-vectorize -fno-tree-slp-vectorize\n";
#else
    out << "vectorization_policy=compiler-specific strict scalar flags\n";
#endif
#else
    out << "strict_scalar=0\n";
    out << "vectorization_policy=compiler default for Release build\n";
#endif
    out << "sve_detected_at_configure=" << AXHEL_BENCHMARK_SVE_DETECTED << '\n';
    out << "axhel_cpu_target=" << AXHEL_BENCHMARK_CPU_TARGET << '\n';
    out << "compiler=" << compiler_name() << '\n';
    out << "cpu_model=" << read_cpu_model() << '\n';
    out << "hostname=" << hostname() << '\n';
    out << "linux_kernel=" << kernel_release() << '\n';
    out << "sve_vector_bits=" << sve_vector_bits() << '\n';
    out << "pinned_cpu=" << pinned_cpu << '\n';
    out << "repetitions=" << cfg.repetitions << '\n';
    out << "warmup=" << cfg.warmup << '\n';
    out << "batch=" << cfg.batch << '\n';
    out << "pmu_profile=" << cfg.pmu_profile << '\n';
    out << "cntvct_enabled=" << (cfg.read_cntvct ? 1 : 0) << '\n';
    out << "seed=0x" << std::hex << cfg.seed << std::dec << '\n';
    out << "scenario_order=" << (cfg.shuffle_scenarios ? "randomized" : "fixed") << '\n';
    out << "shuffle_seed=" << cfg.shuffle_seed << '\n';
    if (cfg.pmu_profile == "none") {
        out << "counter=PERF_TYPE_HARDWARE/PERF_COUNT_HW_CPU_CYCLES\n";
    } else {
        out << "counter=armv8_pmuv3_0 grouped raw PMUv3 events\n";
        out << "cntfrq_el0=" << read_cntfrq_el0() << '\n';
        out << "pmu_cycles_event=0x11\n";
        out << "pmu_instructions_event=0x08\n";
        if (cfg.pmu_profile == "l1") {
            out << "pmu_l1d_cache_rd_event=0x40\n";
            out << "pmu_l1d_cache_refill_event=0x03\n";
        } else if (cfg.pmu_profile == "l2") {
            out << "pmu_l2d_cache_rd_event=0x50\n";
            out << "pmu_l2d_cache_refill_event=0x17\n";
        } else if (cfg.pmu_profile == "stalls") {
            out << "pmu_stall_backend_event=0x24\n";
            out << "pmu_stall_backend_mem_event=0x4005\n";
            out << "pmu_stall_frontend_event=0x23\n";
        }
    }
    out << "counter_scope=user-space only (exclude_kernel=1, exclude_hv=1)\n";
    out << "outlier_filter=Tukey IQR, fences Q1-1.5*IQR and Q3+1.5*IQR\n";
    out << "quantiles=linear interpolation over sorted samples\n";
    out << "ntt_mode=in-place; NTT table construction excluded from measurement\n";
    out << "input_mode=deterministic random, warm-cache steady-state\n";
    out << "batch_cycles_semantics=raw_cycles is total batch PMU count; cycles/statistics are normalized per invocation\n";
}

class CsvWriters {
public:
    CsvWriters(const std::string& output_dir, bool write_raw, std::uint64_t batch, std::string pmu_profile)
        : write_raw_(write_raw), batch_(batch), pmu_profile_(std::move(pmu_profile)), summary_(output_dir + "/summary.csv") {
        if (!summary_) throw std::runtime_error("cannot create summary.csv");
        summary_ << "backend,kernel,poly_degree,qi_bits,modulus,mod_factor,input_mod_factor,output_mod_factor,"
                    "samples_total,samples_kept,outliers_removed,raw_min_cycles,raw_max_cycles,q1_cycles,q3_cycles,iqr_cycles,"
                    "lower_fence_cycles,upper_fence_cycles,min_cycles,max_cycles,mean_cycles,median_cycles,stddev_cycles,cv_percent,"
                    "mean_cycles_per_coeff,median_cycles_per_coeff,perf_scaled\n";

        if (write_raw_) {
            raw_.open(output_dir + "/raw_samples.csv");
            if (!raw_) throw std::runtime_error("cannot create raw_samples.csv");
            raw_ << "backend,kernel,poly_degree,qi_bits,modulus,mod_factor,input_mod_factor,output_mod_factor,"
                    "repetition,batch,pmu_profile,raw_cycles,time_enabled,time_running,cycles,perf_scaled,"
                    "cntvct_ticks,instructions,ipc,l1d_cache_rd,l1d_cache_refill,l2d_cache_rd,l2d_cache_refill,"
                    "stall_backend,stall_backend_mem,stall_frontend\n";
        }
    }

    void write(const KernelMeta& meta, std::uint64_t degree, unsigned bits, std::uint64_t q,
               const std::vector<RawSample>& samples, const Statistics& st) {
        const auto backend = backend_name();
        if (write_raw_) {
            for (const auto& s : samples) {
                raw_ << backend << ',' << meta.name << ',' << degree << ',' << bits << ',' << q << ','
                     << meta.mod_factor << ',' << meta.input_mod_factor << ',' << meta.output_mod_factor << ','
                     << s.repetition << ',' << batch_ << ',' << pmu_profile_ << ','
                     << s.counter.raw_cycles << ',' << s.counter.time_enabled << ','
                     << s.counter.time_running << ',' << std::fixed << std::setprecision(3) << s.counter.cycles << ','
                     << (s.counter.scaled ? 1 : 0) << ','
                     << s.counter.cntvct_ticks << ','
                     << s.counter.instructions << ','
                     << (s.counter.cycles > 0.0 ? s.counter.instructions / s.counter.cycles : 0.0) << ','
                     << s.counter.l1d_cache_rd << ',' << s.counter.l1d_cache_refill << ','
                     << s.counter.l2d_cache_rd << ',' << s.counter.l2d_cache_refill << ','
                     << s.counter.stall_backend << ',' << s.counter.stall_backend_mem << ','
                     << s.counter.stall_frontend << '\n';
            }
        }

        summary_ << backend << ',' << meta.name << ',' << degree << ',' << bits << ',' << q << ','
                 << meta.mod_factor << ',' << meta.input_mod_factor << ',' << meta.output_mod_factor << ','
                 << st.samples_total << ',' << st.samples_kept << ',' << st.outliers_removed << ','
                 << std::fixed << std::setprecision(3)
                 << st.raw_min << ',' << st.raw_max << ',' << st.q1 << ',' << st.q3 << ',' << st.iqr << ','
                 << st.lower_fence << ',' << st.upper_fence << ',' << st.min << ',' << st.max << ','
                 << st.mean << ',' << st.median << ',' << st.stddev << ',' << st.cv_percent << ','
                 << st.mean / static_cast<double>(degree) << ',' << st.median / static_cast<double>(degree) << ','
                 << (st.any_scaled ? 1 : 0) << '\n';
        summary_.flush();
        if (write_raw_) raw_.flush();
    }

private:
    bool write_raw_{true};
    std::uint64_t batch_{1};
    std::string pmu_profile_{"none"};
    std::ofstream raw_;
    std::ofstream summary_;
};

bool selected(const Config& cfg, const std::string& name) {
    return std::find(cfg.kernels.begin(), cfg.kernels.end(), name) != cfg.kernels.end();
}

void report_one(CsvWriters& csv, const KernelMeta& meta, std::uint64_t degree, unsigned bits, std::uint64_t q,
                const std::vector<RawSample>& samples) {
    const auto st = compute_statistics(samples);
    csv.write(meta, degree, bits, q, samples, st);
    std::cout << "  " << std::left << std::setw(20) << meta.name
              << " median=" << std::right << std::fixed << std::setprecision(1) << st.median
              << " cycles  mean=" << st.mean
              << "  kept=" << st.samples_kept << '/' << st.samples_total
              << "  CV=" << std::setprecision(2) << st.cv_percent << "%\n";
}

void run_scenario(const Config& cfg, CsvWriters& csv, std::uint64_t degree, unsigned bits, std::uint64_t q,
                  std::uint64_t root) {
    const auto n = static_cast<std::size_t>(degree);
    std::vector<std::uint64_t> result(n);

    auto run_binary = [&](const std::string& name, auto&& fn) {
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, name));
        auto op1 = random_vector(n, q, 1, rng);
        auto op2 = random_vector(n, q, 1, rng);
        auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
            [] {},
            [&] { fn(result.data(), op1.data(), op2.data()); },
            [&] { g_sink ^= result[(n / 2)]; });
        report_one(csv, {name, 0, 0, 0}, degree, bits, q, samples);
    };

    if (selected(cfg, "add_vv")) {
        run_binary("add_vv", [&](auto* r, const auto* a, const auto* b) { EltwiseAddMod(r, a, b, degree, q); });
    }
    if (selected(cfg, "sub_vv")) {
        run_binary("sub_vv", [&](auto* r, const auto* a, const auto* b) { EltwiseSubMod(r, a, b, degree, q); });
    }

    auto run_scalar = [&](const std::string& name, auto&& fn) {
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, name));
        auto op1 = random_vector(n, q, 1, rng);
        const auto scalar = bounded_random(rng, q);
        auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
            [] {},
            [&] { fn(result.data(), op1.data(), scalar); },
            [&] { g_sink ^= result[(n / 2)]; });
        report_one(csv, {name, 0, 0, 0}, degree, bits, q, samples);
    };

    if (selected(cfg, "add_vs")) {
        run_scalar("add_vs", [&](auto* r, const auto* a, std::uint64_t b) { EltwiseAddMod(r, a, b, degree, q); });
    }
    if (selected(cfg, "sub_vs")) {
        run_scalar("sub_vs", [&](auto* r, const auto* a, std::uint64_t b) { EltwiseSubMod(r, a, b, degree, q); });
    }

    for (const std::uint64_t factor : {1ULL, 2ULL, 4ULL}) {
        const std::string name = "mul_f" + std::to_string(factor);
        if (!selected(cfg, name)) continue;
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, name));
        auto op1 = random_vector(n, q, factor, rng);
        auto op2 = random_vector(n, q, factor, rng);
        auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
            [] {},
            [&] { EltwiseMulMod(result.data(), op1.data(), op2.data(), degree, q, factor); },
            [&] { g_sink ^= result[(n / 2)]; });
        report_one(csv, {name, factor, 0, 0}, degree, bits, q, samples);
    }

    for (const std::uint64_t factor : {1ULL, 2ULL, 4ULL, 8ULL}) {
        const std::string name = "fma_f" + std::to_string(factor);
        if (!selected(cfg, name)) continue;
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, name));
        auto op1 = random_vector(n, q, factor, rng);
        auto op3 = random_vector(n, q, factor, rng);
        const auto scalar = random_factor_value(rng, q, factor);
        auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
            [] {},
            [&] { EltwiseFMAMod(result.data(), op1.data(), scalar, op3.data(), degree, q, factor); },
            [&] { g_sink ^= result[(n / 2)]; });
        report_one(csv, {name, factor, 0, 0}, degree, bits, q, samples);
    }

    struct ReduceCase { const char* name; std::uint64_t in_factor; std::uint64_t out_factor; bool arbitrary; };
    for (const auto rc : {
            ReduceCase{"reduce_2to1", 2, 1, false},
            ReduceCase{"reduce_4to1", 4, 1, false},
            ReduceCase{"reduce_4to2", 4, 2, false},
            ReduceCase{"reduce_barrett", 0, 1, true}}) {
        if (!selected(cfg, rc.name)) continue;
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, rc.name));
        std::vector<std::uint64_t> op(n);
        if (rc.arbitrary) {
            for (auto& x : op) x = rng();
        } else {
            op = random_vector(n, q, rc.in_factor, rng);
        }
        const std::uint64_t public_in_factor = rc.arbitrary ? q : rc.in_factor;
        auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
            [] {},
            [&] { EltwiseReduceMod(result.data(), op.data(), degree, q, public_in_factor, rc.out_factor); },
            [&] { g_sink ^= result[(n / 2)]; });
        report_one(csv, {rc.name, 0, public_in_factor, rc.out_factor}, degree, bits, q, samples);
    }

    const bool need_ntt = selected(cfg, "ntt_forward") || selected(cfg, "ntt_forward_lazy") ||
                          selected(cfg, "ntt_inverse") || selected(cfg, "ntt_inverse_lazy");
    if (need_ntt) {
        NTT ntt(degree, q, root); // table generation explicitly outside measured regions
        std::mt19937_64 rng(scenario_seed(cfg.seed, degree, bits, "ntt_base"));
        auto base = random_vector(n, q, 1, rng);
        auto inverse_base = base;
        ntt.ComputeForward(inverse_base.data(), inverse_base.data(), 1, 1);
        std::vector<std::uint64_t> work(n);

        auto run_ntt = [&](const std::string& name, const std::vector<std::uint64_t>& input,
                           std::uint64_t in_factor, std::uint64_t out_factor, bool forward) {
            if (!selected(cfg, name)) return;
            auto samples = measure_kernel(cfg.warmup, cfg.repetitions, cfg.batch, cfg.pmu_profile, cfg.read_cntvct,
                [&] { std::copy(input.begin(), input.end(), work.begin()); },
                [&] {
                    if (forward) ntt.ComputeForward(work.data(), work.data(), in_factor, out_factor);
                    else ntt.ComputeInverse(work.data(), work.data(), in_factor, out_factor);
                },
                [&] { g_sink ^= work[(n / 2)]; });
            report_one(csv, {name, 0, in_factor, out_factor}, degree, bits, q, samples);
        };

        run_ntt("ntt_forward", base, 1, 1, true);
        run_ntt("ntt_forward_lazy", base, 1, 4, true);
        run_ntt("ntt_inverse", inverse_base, 1, 1, false);
        run_ntt("ntt_inverse_lazy", inverse_base, 1, 2, false);
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config cfg = parse_args(argc, argv);
        const std::string backend = backend_name();
        if (!cfg.expect_backend.empty()) {
            std::string expected = cfg.expect_backend;
            std::transform(expected.begin(), expected.end(), expected.begin(), [](unsigned char c) {
                return static_cast<char>(std::tolower(c));
            });
            std::string actual = backend;
            std::transform(actual.begin(), actual.end(), actual.begin(), [](unsigned char c) {
                return static_cast<char>(std::tolower(c));
            });
            if (expected != actual) {
                throw std::runtime_error("expected backend " + cfg.expect_backend + " but binary uses " + backend);
            }
        }

        const int pinned_cpu = pin_cpu(cfg.cpu);
        ensure_dir(cfg.output_dir);

        // Fail before a long run if the selected PMU mode is unavailable.
        {
            volatile std::uint64_t probe_sink = 0;
            CounterReading reading;
            if (cfg.pmu_profile == "none") {
                PerfCycles probe;
                probe.start();
                for (std::uint64_t i = 0; i < 1000000; ++i) {
                    probe_sink += i ^ (probe_sink >> 1U);
                }
                reading = probe.stop();
            } else {
                PerfGroup probe(cfg.pmu_profile);
                probe.start();
                for (std::uint64_t i = 0; i < 1000000; ++i) {
                    probe_sink += i ^ (probe_sink >> 1U);
                }
                reading = probe.stop();
                if (reading.instructions == 0.0) {
                    throw std::runtime_error("PMU diagnostic instructions counter returned zero");
                }
            }
            g_sink ^= probe_sink;
            if (reading.raw_cycles == 0) {
                throw std::runtime_error("PMU cycle counter returned zero during busy-loop probe");
            }
        }

        write_metadata(cfg, pinned_cpu, cfg.output_dir);
        CsvWriters csv(cfg.output_dir, cfg.write_raw, cfg.batch, cfg.pmu_profile);

        std::cout << "AXHEL benchmark backend: " << backend << '\n'
                  << "CPU target: " << AXHEL_BENCHMARK_CPU_TARGET << '\n'
                  << "Pinned CPU: " << pinned_cpu << '\n'
                  << "Repetitions: " << cfg.repetitions << ", warm-up: " << cfg.warmup
                  << ", batch: " << cfg.batch << '\n'
                  << "PMU profile: " << cfg.pmu_profile
                  << ", CNTVCT: " << (cfg.read_cntvct ? "on" : "off") << '\n'
                  << "Output: " << cfg.output_dir << "\n";

        struct Scenario {
            std::uint64_t degree;
            unsigned bits;
        };

        std::vector<Scenario> scenarios;
        scenarios.reserve(cfg.degrees.size() * cfg.bits.size());
        for (const auto degree : cfg.degrees) {
            for (const auto bits : cfg.bits) {
                scenarios.push_back({degree, bits});
            }
        }

        if (cfg.shuffle_scenarios) {
            std::mt19937_64 shuffle_rng(cfg.shuffle_seed);
            std::shuffle(scenarios.begin(), scenarios.end(), shuffle_rng);
        }

        std::size_t scenario_index = 0;
        for (const auto& scenario : scenarios) {
            const auto degree = scenario.degree;
            const auto bits = scenario.bits;
            const auto q = find_ntt_prime(degree, bits);
            const auto root = find_primitive_2n_root(degree, q);
            std::cout << "\nScenario " << (++scenario_index) << '/' << scenarios.size()
                      << " N=" << degree << ", qi_bits=" << bits << ", q=" << q << '\n';
            run_scenario(cfg, csv, degree, bits, q, root);
        }

        std::cout << "\nCompleted. sink=" << g_sink << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "AXHEL benchmark error: " << e.what() << '\n';
        return 1;
    }
}
