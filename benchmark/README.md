# AXHEL Benchmark

AXHEL includes a benchmarking framework for evaluating the performance of its modular arithmetic and Number Theoretic Transform (NTT) kernels.

The benchmark exercises the AXHEL public API and compares three separately compiled configurations under equivalent build and execution conditions: ScalarStrict, Native, and SVE.

Measurements are based on Linux hardware performance counters and use deterministic benchmark inputs and parameter generation.

## Benchmark configurations

The benchmark compares three AXHEL configurations:

| Configuration | Description |
|---|---|
| `ScalarStrict` | Native AXHEL implementation with compiler auto-vectorization disabled |
| `Native` | Native AXHEL implementation using the normal Release compiler optimization policy |
| `SVE` | Explicit AXHEL SVE implementation selected through normal SVE detection |

All configurations use `CMAKE_BUILD_TYPE=Release` and the same `AXHEL_CPU` target.

### ScalarStrict

`ScalarStrict` provides a scalar reference while preserving the normal Release optimization level and CPU-specific tuning.

The AXHEL Native implementation is selected, explicit SVE source files are excluded, and compiler-generated loop and SLP vectorization are disabled for both the AXHEL library and benchmark executable.

For GCC, the following options are used:

```text
-fno-tree-vectorize
-fno-tree-loop-vectorize
-fno-tree-slp-vectorize
```

For Clang:

```text
-fno-vectorize
-fno-slp-vectorize
```

Other compiler optimizations, including instruction scheduling, inlining, loop unrolling, and CPU-specific tuning, remain enabled.

The benchmark runner verifies the generated compilation commands before execution to ensure that the vectorization-disable flags are applied to the AXHEL library sources and that explicit SVE translation units are not included.

### Native

`Native` uses the AXHEL Native backend with the normal Release compiler optimization policy.

The compiler is therefore free to auto-vectorize Native source code when appropriate.

### SVE

`SVE` uses normal AXHEL SVE capability detection and the public AXHEL dispatcher selects the explicit SVE implementations.

The complete three-configuration comparison therefore requires an SVE-capable compiler and target.

All three configurations invoke the same public AXHEL API. Internal implementation functions are not called directly by the benchmark.

## CMake configuration options

The following CMake options control benchmark builds:

| Option | Default | Description |
|---|:---:|---|
| `AXHEL_BENCHMARK` | `OFF` | Enables the AXHEL benchmark target |
| `AXHEL_BENCHMARK_FORCE_NATIVE` | `OFF` | Forces the Native backend for a benchmark build |
| `AXHEL_BENCHMARK_FORCE_SCALAR` | `OFF` | Forces the strict scalar benchmark configuration |
| `AXHEL_CPU` | `native` | Selects the Arm CPU target passed to the compiler through `-mcpu` |

`AXHEL_BENCHMARK_FORCE_NATIVE` and `AXHEL_BENCHMARK_FORCE_SCALAR` are mutually exclusive and are valid only when `AXHEL_BENCHMARK=ON`.

Normal AXHEL builds are unaffected because benchmarking is disabled by default.

## Benchmarked kernels

The benchmark covers the public modular arithmetic and NTT APIs.

| Kernel identifier | AXHEL operation |
|---|---|
| `add_vv` | `EltwiseAddMod`, vector-vector |
| `add_vs` | `EltwiseAddMod`, vector-scalar |
| `sub_vv` | `EltwiseSubMod`, vector-vector |
| `sub_vs` | `EltwiseSubMod`, vector-scalar |
| `mul_f1` | `EltwiseMulMod`, input modulus factor 1 |
| `mul_f2` | `EltwiseMulMod`, input modulus factor 2 |
| `mul_f4` | `EltwiseMulMod`, input modulus factor 4 |
| `fma_f1` | `EltwiseFMAMod`, input modulus factor 1 |
| `fma_f2` | `EltwiseFMAMod`, input modulus factor 2 |
| `fma_f4` | `EltwiseFMAMod`, input modulus factor 4 |
| `fma_f8` | `EltwiseFMAMod`, input modulus factor 8 |
| `reduce_2to1` | `EltwiseReduceMod`, factor 2 to factor 1 |
| `reduce_4to1` | `EltwiseReduceMod`, factor 4 to factor 1 |
| `reduce_4to2` | `EltwiseReduceMod`, factor 4 to factor 2 |
| `reduce_barrett` | `EltwiseReduceMod`, arbitrary-uint64 Barrett path |
| `ntt_forward` | `NTT::ComputeForward`, normalized output |
| `ntt_forward_lazy` | `NTT::ComputeForward`, lazy output |
| `ntt_inverse` | `NTT::ComputeInverse`, normalized output |
| `ntt_inverse_lazy` | `NTT::ComputeInverse`, lazy output |

## Default benchmark configuration

Running:

```bash
./benchmark/scripts/run-benchmark.sh
```

uses the following defaults:

| Setting | Default |
|---|---|
| Build type | `Release` |
| C++ compiler | `$CXX`, or `g++` if unset |
| AXHEL CPU target | `native` |
| Measured repetitions | `100` |
| Warm-up repetitions | `10` |
| Batch size | `1` |
| PMU profile | `none` |
| CNTVCT measurement | disabled |
| CPU affinity | `auto` |
| Polynomial degrees | `2048,4096,8192,16384,32768` |
| Modulus bit widths | `20,25,30,40,50,52,54,56,58,60,61,62` |
| Kernels | `all` |
| Scenario order | fixed |
| Raw samples | enabled |
| Output directory | `benchmark/results/<timestamp>/` |

With `--cpu auto`, the benchmark pins itself to the first CPU available in the current process affinity mask.

## Benchmark parameter matrix

For each selected polynomial degree and modulus bit width, the benchmark deterministically searches for an NTT-compatible prime modulus satisfying:

```text
q = 1 mod 2N
```

with exactly the requested bit width.

The same modulus is used by the arithmetic and NTT kernels associated with that scenario.

Every measured kernel invocation operates on `N` coefficients.

Input data are generated deterministically so that corresponding ScalarStrict, Native, and SVE measurements use equivalent datasets.

## Measurement methodology

CPU cycles are measured directly through the Linux `perf_event_open()` interface using:

```text
PERF_TYPE_HARDWARE
PERF_COUNT_HW_CPU_CYCLES
```

The hardware counter measures user-space execution only; kernel and hypervisor execution are excluded.

The benchmark records both the raw counter value and the Linux `time_enabled` and `time_running` values. If an event is multiplexed, the cycle count is automatically scaled and the corresponding output is marked through the `perf_scaled` field.

Warm-up invocations are executed before measured repetitions.

NTT table construction, input preparation, and input restoration are performed outside the measured interval. NTT operations are measured in-place.

## Statistical processing

Measured cycle counts are processed independently for every backend, kernel, polynomial degree, and modulus configuration.

Outliers are identified using Tukey's IQR method:

```text
lower fence = Q1 - 1.5 * IQR
upper fence = Q3 + 1.5 * IQR
```

Samples outside the fences are excluded from the final statistics.

The generated summaries include the number of retained and rejected samples, raw and filtered ranges, Q1, Q3, IQR, mean, median, sample standard deviation, coefficient of variation, and cycles per coefficient.

## Benchmark scripts

Benchmark utilities are located in [`benchmark/scripts/`](scripts/).

| Script | Description |
|---|---|
| `run-benchmark.sh` | Main benchmark runner; configures, builds, executes, and compares ScalarStrict, Native, and SVE |
| `compare_results.py` | Merges the three backend summaries and computes speedups and cycle reductions |
| `run-all-kernels-screening.sh` | Runs repeated randomized all-kernel benchmark campaigns |
| `run-all-kernels-definitive.sh` | Runs repeated all-kernel campaigns with separate element-wise and NTT measurement settings |

`run-benchmark.sh` is the main interface for general benchmarking and custom benchmark configurations.

The campaign scripts are higher-level wrappers around `run-benchmark.sh` and define their own repetition counts, batching policies, CPU affinity, parameter matrices, and scenario ordering.

## Run the benchmark

From the AXHEL repository root:

```bash
chmod +x benchmark/scripts/*.sh benchmark/scripts/*.py
./benchmark/scripts/run-benchmark.sh
```

The runner performs clean builds and benchmark executions in the following order:

```text
ScalarStrict
Native
SVE
```

It uses separate build directories:

```text
build-benchmark-scalar-strict/
build-benchmark-native/
build-benchmark-sve/
```

After the three executions, `compare_results.py` automatically generates the combined speedup report.

## Runner options

The default configuration can be changed through the following options:

| Option | Description |
|---|---|
| `--repetitions N` | Number of measured repetitions per scenario |
| `--warmup N` | Number of warm-up repetitions |
| `--batch N` | Number of kernel invocations included in each PMU sample |
| `--pmu-profile P` | PMU profile: `none`, `l1`, `l2`, or `stalls` |
| `--cntvct` | Enables `CNTVCT_EL0` measurement for diagnostic PMU profiles |
| `--cpu auto\|off\|ID` | Controls CPU affinity |
| `--degrees LIST` | Comma-separated polynomial degrees |
| `--bits LIST` | Comma-separated modulus bit widths |
| `--kernels LIST\|all` | Selects the kernels to benchmark |
| `--axhel-cpu CPU_TARGET` | Sets the value passed to `AXHEL_CPU` |
| `--shuffle-scenarios` | Randomizes the scenario execution order |
| `--shuffle-seed N` | Sets the deterministic scenario-order seed |
| `--output DIR` | Selects the benchmark output directory |
| `--no-raw` | Disables generation of `raw_samples.csv` |

For example:

```bash
./benchmark/scripts/run-benchmark.sh \
    --cpu 10 \
    --repetitions 1000 \
    --warmup 100
```

A subset of the benchmark matrix can be selected with:

```bash
./benchmark/scripts/run-benchmark.sh \
    --degrees 8192,16384,32768 \
    --bits 40,50,52,60,61,62
```

Specific kernels can be selected with:

```bash
./benchmark/scripts/run-benchmark.sh \
    --kernels add_vv,sub_vv,mul_f1,fma_f1,reduce_barrett,ntt_forward,ntt_inverse
```

## CPU affinity

CPU affinity is controlled through:

```text
--cpu auto
--cpu off
--cpu ID
```

`auto`, the default, pins the benchmark to the first CPU allowed by the current process affinity mask.

`off` disables benchmark-controlled affinity.

A numeric CPU identifier pins the benchmark to that specific CPU, provided that it is available in the current affinity mask.

For reproducible measurements, a dedicated or isolated CPU can be selected explicitly.

## Batch measurements

The `--batch` option controls the number of consecutive kernel invocations included in one PMU measurement.

The reported cycle statistics are normalized per invocation.

The default is:

```text
--batch 1
```

Batch values greater than one are supported only for element-wise kernels.

NTT operations are in-place and therefore require:

```text
--batch 1
```

## Scenario ordering

By default, `(poly_degree, qi_bits)` scenarios are executed in a fixed order.

The scenario order can be randomized with:

```bash
./benchmark/scripts/run-benchmark.sh \
    --shuffle-scenarios \
    --shuffle-seed 12345
```

Randomization affects only the execution order. The benchmark input data remain deterministic.

Using the same shuffle seed reproduces the same scenario order.

## PMU diagnostic profiles

The default benchmark mode measures CPU cycles only:

```text
--pmu-profile none
```

Additional diagnostic profiles are available:

| Profile | Additional information |
|---|---|
| `l1` | Instructions and L1 data-cache access/refill events |
| `l2` | Instructions and L2 data-cache access/refill events |
| `stalls` | Instructions and backend/frontend stall events |

Diagnostic profiles use grouped Arm PMU events and require:

```text
--batch 1
```

They expect the Linux Arm PMU interface exposed as `armv8_pmuv3_0`. Availability of individual PMU events may depend on the target processor.

On AArch64, `CNTVCT_EL0` ticks can additionally be recorded with:

```bash
./benchmark/scripts/run-benchmark.sh \
    --pmu-profile l1 \
    --cntvct
```

`--cntvct` is available only together with a diagnostic PMU profile.

## Environment variables

The benchmark runner recognizes the following environment variables:

| Variable | Default | Description |
|---|:---:|---|
| `CXX` | `g++` | C++ compiler used for the three benchmark builds |
| `AXHEL_BENCH_JOBS` | CMake default | Number of parallel build jobs |

For example:

```bash
CXX=g++ AXHEL_BENCH_JOBS=8 \
./benchmark/scripts/run-benchmark.sh \
    --axhel-cpu neoverse-v2
```

## Output files

A standard benchmark execution creates:

```text
benchmark/results/<timestamp>/
├── scalar-strict/
│   ├── metadata.txt
│   ├── raw_samples.csv
│   └── summary.csv
├── native/
│   ├── metadata.txt
│   ├── raw_samples.csv
│   └── summary.csv
├── sve/
│   ├── metadata.txt
│   ├── raw_samples.csv
│   └── summary.csv
├── speedup.csv
└── run-config.txt
```

`metadata.txt` records information about the backend, compiler, CPU target, hardware, affinity, benchmark configuration, PMU configuration, SVE vector length, and statistical processing.

`raw_samples.csv` contains the individual PMU measurements. It is omitted when `--no-raw` is specified.

`summary.csv` contains the IQR-filtered statistics for every benchmark scenario.

`speedup.csv` combines the ScalarStrict, Native, and SVE summaries and reports the performance comparisons between the three configurations.

`run-config.txt` records the configuration used by the benchmark runner.

## Speedup comparisons

The combined `speedup.csv` reports three comparisons:

```text
ScalarStrict / SVE
ScalarStrict / Native
Native / SVE
```

`ScalarStrict / SVE` measures the improvement of the explicit SVE implementation relative to the same AXHEL public implementation compiled without compiler-generated vectorization.

`ScalarStrict / Native` shows the contribution of the compiler's normal Release auto-vectorization policy to the Native implementation.

`Native / SVE` compares the explicit AXHEL SVE implementation with the normally optimized Native implementation.

Each comparison is reported using both IQR-filtered median and mean cycle counts, together with the corresponding cycle-reduction percentages.

## Campaign scripts

Two additional scripts provide repeated benchmark campaigns:

```text
benchmark/scripts/run-all-kernels-screening.sh
benchmark/scripts/run-all-kernels-definitive.sh
```

Both scripts repeatedly invoke `run-benchmark.sh`, randomize the benchmark scenario order using deterministic shuffle seeds, and store each execution in a separate result directory.

They also record the Git commit, working-tree status, source diff, campaign configuration, shuffle seeds, and execution logs.

The screening script runs a single all-kernel configuration for each campaign iteration.

The definitive script separates element-wise/reduction kernels from NTT kernels so that different batch sizes can be used while preserving the in-place NTT constraint.

These scripts contain predefined campaign parameters, including CPU affinity and benchmark matrices. They should therefore be reviewed and, when necessary, adapted to the target system before use.

For general benchmarking and custom configurations, use `run-benchmark.sh`.

## PMU permissions

Hardware performance-counter access depends on the Linux system configuration.

If PMU access is denied, inspect:

```bash
cat /proc/sys/kernel/perf_event_paranoid
```

and configure the system according to the applicable local security policy.

The benchmark performs a PMU probe before starting the complete measurement campaign and terminates early if the selected counters are unavailable.

## Normal AXHEL builds

Benchmarking is disabled by default:

```text
AXHEL_BENCHMARK=OFF
```

When benchmarking is disabled, the benchmark executable and benchmark-specific backend overrides are not included in the normal AXHEL build.