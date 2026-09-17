# AXHEL public-kernel benchmark: ScalarStrict vs Native vs SVE

This benchmark measures Linux hardware PMU CPU cycles around the AXHEL **public API**
and compares three separately compiled Release builds of exactly the same benchmark
workload.

## The three baselines

### 1. ScalarStrict

`ScalarStrict` is the strict scalar reference used to quantify the benefit of
vectorization itself.

- normal SVE capability detection is still executed at CMake configure time;
- the effective AXHEL backend is forced to the public Native path;
- explicit AXHEL SVE translation units are not compiled;
- the compiler auto-vectorizer is disabled on the **AXHEL library sources** and,
  as an additional guard, on the benchmark driver itself;
- `-O3`, the selected `-mcpu=${AXHEL_CPU}` target, scalar instruction scheduling,
  inlining, loop unrolling and other non-vector optimizations remain enabled.

For GCC the build appends:

    -fno-tree-vectorize
    -fno-tree-loop-vectorize
    -fno-tree-slp-vectorize

For Clang it appends:

    -fno-vectorize
    -fno-slp-vectorize

The runner checks `compile_commands.json` before execution to verify that the
vectorization-disable flags are actually attached to AXHEL source files and that
no explicit SVE source file is present in the ScalarStrict build.

This is intentionally **not** implemented by changing the CPU target to an older
architecture without SIMD. All three builds remain tuned for the same processor,
so the comparison does not mix the effect of vectorization with a different
instruction-set/tuning target.

### 2. Native

`Native` uses the AXHEL public Native backend but leaves the compiler's normal
Release optimization policy unchanged. Therefore GCC/Clang may auto-vectorize
scalar source code when it considers this profitable.

### 3. SVE

`SVE` uses normal AXHEL automatic SVE detection and the public dispatcher selects
the explicit AXHEL SVE kernels.

All three executables call the same public API surface. Private implementation
symbols are never called directly by the benchmark.

## Interpretation of the comparisons

The combined `speedup.csv` reports three complementary comparisons.

Primary vectorization comparison:

    speedup_median_scalar_over_sve = median_scalar_strict / median_sve

This answers: **how much faster is the explicit SVE implementation than the same
AXHEL public kernel compiled without compiler vectorization?**

Compiler auto-vectorization contribution:

    speedup_median_scalar_over_native = median_scalar_strict / median_native

This shows how much the ordinary compiler optimization policy alone improves the
Native implementation.

Explicit SVE versus normal Native:

    speedup_median_native_over_sve = median_native / median_sve

This answers: **how much does AXHEL's explicit SVE implementation improve on the
best normal Native Release build?**

The same three comparisons are also emitted using the IQR-filtered mean, together
with the corresponding cycle-reduction percentages.

## Default HE-representative matrix

Polynomial degree:

    2048, 4096, 8192, 16384, 32768

Coefficient-modulus width:

    20, 25, 30, 40, 50, 52, 54, 56, 58, 60, 61, 62 bits

For every `(poly_degree, qi_bits)` pair the benchmark deterministically searches
for a prime with exactly the requested bit width satisfying
`q = 1 (mod 2*poly_degree)`. The same NTT-friendly `q` is used by all kernels in
that scenario, so arithmetic and NTT measurements refer to a consistent HE
parameter point.

Every measured kernel invocation processes exactly `N = poly_degree`
coefficients. The default execution uses 10 warm-up invocations followed by 100
measured repetitions per kernel/scenario.

## Public kernels covered

- `EltwiseAddMod`: vector/vector and vector/scalar
- `EltwiseSubMod`: vector/vector and vector/scalar
- `EltwiseMulMod`: input modulus factors 1, 2, 4
- `EltwiseFMAMod`: input modulus factors 1, 2, 4, 8
- `EltwiseReduceMod`: 2->1, 4->1, 4->2 and arbitrary-uint64 Barrett path
- `NTT::ComputeForward`: normalized and lazy output (`output_mod_factor=4`)
- `NTT::ComputeInverse`: normalized and lazy output (`output_mod_factor=2`)

NTT table construction and input restoration are outside the measured interval.
NTT execution is in-place to measure the transform rather than an out-of-place
copy. Inputs are deterministic and warm-ups establish a steady-state/warm-cache
microbenchmark regime.

## Cycle measurement

The executable uses Linux `perf_event_open()` directly with:

    type   = PERF_TYPE_HARDWARE
    config = PERF_COUNT_HW_CPU_CYCLES

Kernel/hypervisor execution is excluded. The benchmark records raw event values
plus `time_enabled/time_running`; if Linux multiplexes the event, the cycle value
is scaled and the CSV marks `perf_scaled=1`.

By default the process pins itself to the first CPU allowed by its affinity mask.
Use `--cpu ID` to choose a specific core or `--cpu off` to disable pinning.

If PMU access is denied, inspect:

    cat /proc/sys/kernel/perf_event_paranoid

and configure the machine according to local system policy.

## IQR outlier filtering and statistics

For each independent `(backend, kernel, poly_degree, qi_bits)` sample set:

1. sort measured cycle counts;
2. compute Q1 and Q3 at the 25th and 75th percentiles using linear interpolation;
3. compute `IQR = Q3 - Q1`;
4. remove samples outside `[Q1 - 1.5*IQR, Q3 + 1.5*IQR]`;
5. compute final statistics on retained samples.

Each summary contains raw min/max, filtered min/max, mean, median, sample standard
deviation, coefficient of variation, Q1/Q3/IQR and cycles per coefficient.

## Build and run all three backends

From the repository root:

    chmod +x benchmark/scripts/*.sh benchmark/scripts/*.py
    ./benchmark/scripts/run-benchmark.sh

The runner performs, in order:

    ScalarStrict clean build -> benchmark
    Native       clean build -> benchmark
    SVE          clean build -> benchmark
    merge + speedup computation

Results are created under:

    benchmark/results/YYYYMMDD-HHMMSS/
        scalar-strict/metadata.txt
        scalar-strict/raw_samples.csv
        scalar-strict/summary.csv
        native/metadata.txt
        native/raw_samples.csv
        native/summary.csv
        sve/metadata.txt
        sve/raw_samples.csv
        sve/summary.csv
        speedup.csv
        run-config.txt

A publication-oriented run pinned to a chosen isolated core can be launched as:

    ./benchmark/scripts/run-benchmark.sh \
        --cpu 10 \
        --repetitions 1000 \
        --warmup 100

To benchmark only the main normalized kernels:

    ./benchmark/scripts/run-benchmark.sh \
        --kernels add_vv,sub_vv,mul_f1,fma_f1,reduce_barrett,ntt_forward,ntt_inverse

You can also override the parameter matrix, for example:

    ./benchmark/scripts/run-benchmark.sh \
        --degrees 8192,16384,32768 \
        --bits 40,50,52,60,61,62

## Recommended paper terminology

A concise description is:

> We compare three Release configurations: a strict scalar baseline with compiler
> loop and SLP vectorization disabled, the normally optimized native backend, and
> the explicit SVE backend. All configurations target the same processor and are
> evaluated through the same public AXHEL API.

This terminology avoids claiming that the physical processor lacks SIMD/SVE. The
hardware is identical across runs; what is removed in `ScalarStrict` is all
compiler-generated vectorization together with AXHEL's explicit SVE backend.
