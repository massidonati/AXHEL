#pragma once
#include "reference.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif

namespace axhel_test {

struct ModulusCase {
    unsigned bits;
    std::uint64_t q;
    unsigned variant; // 0 = primary/high, 1 = low-quarter, 2 = mid-range
};

inline unsigned bit_width(std::uint64_t x) {
    unsigned b = 0;
    while (x) { ++b; x >>= 1U; }
    return b;
}

inline bool is_prime64(std::uint64_t n) {
    if (n < 2) return false;
    for (std::uint64_t p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL}) {
        if (n % p == 0) return n == p;
    }
    std::uint64_t d = n - 1;
    unsigned s = 0;
    while ((d & 1U) == 0) { d >>= 1U; ++s; }

    // Deterministic Miller-Rabin bases for uint64_t.
    for (std::uint64_t a : {2ULL,325ULL,9375ULL,28178ULL,450775ULL,9780504ULL,1795265022ULL}) {
        if (a % n == 0) continue;
        std::uint64_t x = pow_mod(a % n, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;
        for (unsigned r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) { witness = false; break; }
        }
        if (witness) return false;
    }
    return true;
}

inline std::uint64_t align_ntt_candidate_down(std::uint64_t x, std::uint64_t step) {
    // Return the largest q <= x such that q == 1 (mod step).
    return x - ((x - 1) % step);
}

inline std::uint64_t find_ntt_prime(std::uint64_t degree, unsigned bits) {
    if (bits < 2 || bits >= 63) throw std::invalid_argument("test helper supports 2 <= bits < 63");
    const std::uint64_t step = 2 * degree;
    const std::uint64_t hi = (std::uint64_t{1} << bits) - 1;
    const std::uint64_t lo = (std::uint64_t{1} << (bits - 1));
    std::uint64_t q = align_ntt_candidate_down(hi, step);
    while (q >= lo && q > step) {
        if (is_prime64(q)) return q;
        if (q < step) break;
        q -= step;
    }
    throw std::runtime_error("unable to find NTT-friendly prime");
}

// Deterministically select an NTT-friendly prime near a requested position in
// the bit-width interval [2^(bits-1), 2^bits).  No random state is involved.
inline std::uint64_t find_ntt_prime_near(std::uint64_t degree, unsigned bits,
                                         std::uint64_t numerator,
                                         std::uint64_t denominator) {
    if (bits < 2 || bits >= 63 || denominator == 0 || numerator > denominator)
        throw std::invalid_argument("invalid deterministic NTT-prime request");

    const std::uint64_t step = 2 * degree;
    const std::uint64_t lo = (std::uint64_t{1} << (bits - 1));
    const std::uint64_t hi = (std::uint64_t{1} << bits) - 1;
    const u128 span = static_cast<u128>(hi) - lo;
    const std::uint64_t target = lo + static_cast<std::uint64_t>(span * numerator / denominator);

    // Search downward first so the result is stable across platforms.
    std::uint64_t q = align_ntt_candidate_down(target, step);
    while (q >= lo && q > step) {
        if (is_prime64(q)) return q;
        if (q < lo + step) break;
        q -= step;
    }

    // Deterministic fallback upward from the aligned target.
    q = align_ntt_candidate_down(target, step);
    if (q < target) {
        if (q > std::numeric_limits<std::uint64_t>::max() - step)
            throw std::runtime_error("NTT-prime search overflow");
        q += step;
    }
    while (q <= hi && bit_width(q) == bits) {
        if (is_prime64(q)) return q;
        if (q > hi - step) break;
        q += step;
    }

    throw std::runtime_error("unable to find NTT-friendly prime near requested position");
}

inline std::uint64_t find_primitive_2n_root(std::uint64_t degree, std::uint64_t q) {
    const std::uint64_t order = 2 * degree;
    const std::uint64_t exponent = (q - 1) / order;
    for (std::uint64_t a = 2; a < 100000; ++a) {
        const auto psi = pow_mod(a, exponent, q);
        if (psi != 1 &&
            pow_mod(psi, degree, q) == q - 1 &&
            pow_mod(psi, order, q) == 1) {
            return psi;
        }
    }
    throw std::runtime_error("unable to find primitive root");
}

inline std::uint64_t checked_bound(std::uint64_t q, std::uint64_t factor) {
    const u128 x = static_cast<u128>(q) * factor;
    if (x > std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error("q*factor exceeds uint64_t");
    return static_cast<std::uint64_t>(x);
}

inline std::uint64_t effective_uint64_bound(std::uint64_t q, std::uint64_t factor) {
    const u128 mathematical_bound = static_cast<u128>(q) * factor;
    const u128 domain_size = static_cast<u128>(std::numeric_limits<std::uint64_t>::max()) + 1;
    // random_vector takes an exclusive uint64 bound. UINT64_MAX therefore
    // represents [0, UINT64_MAX), while exact UINT64_MAX is covered explicitly
    // by boundary tests.
    if (mathematical_bound >= domain_size)
        return std::numeric_limits<std::uint64_t>::max();
    return static_cast<std::uint64_t>(mathematical_bound);
}

inline constexpr std::uint64_t test_seed() noexcept {
    // V5 regression is intentionally fixed: every execution uses the exact
    // same pseudo-random streams and parameter matrix.
    return 0x415848454cULL; // "AXHEL"
}

inline std::size_t vector_lanes_u64() {
#ifdef __ARM_FEATURE_SVE
    return static_cast<std::size_t>(svcntd());
#else
    return 2; // useful surrogate for forced-Native test builds
#endif
}

// Single deterministic V5 regression matrix. It is intentionally independent
// of environment variables and therefore identical on every execution.
inline std::vector<std::size_t> element_sizes() {
    const std::size_t L = vector_lanes_u64();
    std::vector<std::size_t> v = {
        0,1,2,3,4,5,6,7,
        31,32,33,
        127,128,129,
        255,256,257,
        1023,1024,1025,
        8191,8192,8193,
        16383,16384,16385,
        32767,32768,32769
    };

    // Explicit SVE-tail/unroll boundaries.  These remain useful in Native
    // builds as ordinary short-vector boundary cases.
    for (std::size_t k : {1ULL,2ULL,3ULL,4ULL,7ULL,8ULL,16ULL}) {
        const std::size_t x = k * L;
        if (x > 0) v.push_back(x - 1);
        v.push_back(x);
        v.push_back(x + 1);
    }

    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

inline std::vector<std::size_t> canary_tail_sizes() {
    const std::size_t L = vector_lanes_u64();
    std::vector<std::size_t> v{0,1,2,3,31,32,33,127,128,129,1023,1024,1025,8191,8192,8193};
    for (std::size_t k : {1ULL,2ULL,3ULL,4ULL,7ULL,8ULL,16ULL}) {
        const std::size_t x = k * L;
        if (x > 0) v.push_back(x - 1);
        v.push_back(x);
        v.push_back(x + 1);
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

// General arithmetic/NTT classes: previous exhaustive set + explicit 62-bit.
inline const std::vector<unsigned>& modulus_bits() {
    static const std::vector<unsigned> v{20,30,40,50,52,55,59,60,61,62};
    return v;
}

// Exercise every public Barrett-dispatch shift used by HE-sized moduli:
// BitWidth(q)=20..62 => Shift=18..60.
inline const std::vector<unsigned>& barrett_modulus_bits() {
    static const std::vector<unsigned> v = [] {
        std::vector<unsigned> out;
        for (unsigned bits = 20; bits <= 62; ++bits) out.push_back(bits);
        return out;
    }();
    return v;
}

inline bool is_critical_barrett_bits(unsigned bits) {
    return bits==50 || bits==52 || bits==59 || bits==60 || bits==61 || bits==62;
}

// One deterministic high-end q for every bit-width 20..62.  For critical HE
// widths, add two distinct q values from lower/middle portions of the interval.
inline const std::vector<ModulusCase>& barrett_modulus_cases() {
    static const std::vector<ModulusCase> cases = [] {
        constexpr std::uint64_t degree = 1024;
        std::vector<ModulusCase> out;
        for (auto bits : barrett_modulus_bits()) {
            const auto high = find_ntt_prime(degree, bits);
            out.push_back({bits, high, 0});
            if (is_critical_barrett_bits(bits)) {
                const auto low = find_ntt_prime_near(degree, bits, 1, 4);
                const auto mid = find_ntt_prime_near(degree, bits, 1, 2);
                if (low != high) out.push_back({bits, low, 1});
                if (mid != high && mid != low) out.push_back({bits, mid, 2});
            }
        }
        return out;
    }();
    return cases;
}

inline int random_repetitions() { return 256; }
inline int ntt_random_repetitions() { return 32; }
inline int polynomial_oracle_repetitions() { return 64; }
inline int shift_random_repetitions() { return 1024; }

inline const std::vector<std::size_t>& ntt_degrees() {
    static const std::vector<std::size_t> v{
        8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768
    };
    return v;
}

inline std::uint64_t bounded_random(std::mt19937_64& rng, std::uint64_t bound) {
    if (bound == 0) throw std::invalid_argument("bounded_random bound must be non-zero");
    // Rejection sampling avoids implementation-defined behavior of
    // std::uniform_int_distribution, making the generated sequence depend only
    // on the standardized mt19937_64 engine and the fixed seed.
    const std::uint64_t threshold = static_cast<std::uint64_t>(-bound) % bound;
    for (;;) {
        const std::uint64_t x = rng();
        if (x >= threshold) return x % bound;
    }
}

inline std::vector<std::uint64_t> random_vector(std::size_t n, std::uint64_t bound,
                                                std::mt19937_64& rng) {
    std::vector<std::uint64_t> out(n);
    for (auto& x : out) x = bounded_random(rng, bound);
    return out;
}

inline std::vector<std::uint64_t> boundary_values(std::uint64_t q, std::uint64_t factor) {
    const auto B = checked_bound(q, factor);
    std::vector<std::uint64_t> v{0,1};
    auto add_if = [&](u128 x) {
        if (x < B && x <= std::numeric_limits<std::uint64_t>::max())
            v.push_back(static_cast<std::uint64_t>(x));
    };

    for (std::uint64_t k = 1; k < factor; ++k) {
        const u128 t = static_cast<u128>(k) * q;
        if (t >= 2) add_if(t - 2);
        if (t >= 1) add_if(t - 1);
        add_if(t);
        add_if(t + 1);
    }

    if (B >= 2) add_if(static_cast<u128>(B)-2);
    if (B >= 1) add_if(static_cast<u128>(B)-1);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

inline std::vector<std::uint64_t> boundary_values_uint64_domain(
    std::uint64_t q, std::uint64_t factor)
{
    const u128 mathematical_bound = static_cast<u128>(q) * factor;
    const u128 uint64_limit = static_cast<u128>(std::numeric_limits<std::uint64_t>::max()) + 1;
    const u128 exclusive_bound = std::min(mathematical_bound, uint64_limit);

    std::vector<std::uint64_t> v{0,1};
    auto add_if = [&](u128 x) {
        if (x < exclusive_bound && x <= std::numeric_limits<std::uint64_t>::max())
            v.push_back(static_cast<std::uint64_t>(x));
    };

    if (q > 1) add_if(q-2);
    add_if(q-1); add_if(q); add_if(static_cast<u128>(q)+1);
    add_if(static_cast<u128>(2)*q-1); add_if(static_cast<u128>(2)*q);
    add_if(static_cast<u128>(2)*q+1);
    add_if(static_cast<u128>(3)*q-1); add_if(static_cast<u128>(3)*q);
    add_if(static_cast<u128>(3)*q+1);
    add_if(static_cast<u128>(4)*q-1); add_if(static_cast<u128>(4)*q);
    add_if(static_cast<u128>(4)*q+1);
    add_if(static_cast<u128>(7)*q-1); add_if(static_cast<u128>(7)*q);
    add_if(static_cast<u128>(7)*q+1);

    if (exclusive_bound > 1) add_if(exclusive_bound-2);
    if (exclusive_bound > 0) add_if(exclusive_bound-1);

    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

inline std::vector<std::uint64_t> make_pattern(std::size_t n,
                                               const std::vector<std::uint64_t>& vals) {
    std::vector<std::uint64_t> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = vals[i % vals.size()];
    return v;
}

} // namespace axhel_test
