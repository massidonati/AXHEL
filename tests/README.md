# AXHEL deterministic extended regression suite — V5

V5 intentionally provides **one test mode only: regression**. There is no smoke,
exhaustive, or stress mode. The regression matrix is fixed and deterministic so
that every execution exercises the same cases in the same order.

## Design goals

- catch functional regressions in both Native and SVE backends;
- test public APIs against independent mathematical oracles;
- directly test critical Native/SVE helpers and NTT butterflies;
- exercise the complete HE-relevant Barrett bit-width range through 62 bits;
- cover every SVE `ShiftRight128LowPart<Shift>` specialization from 1 through 60;
- use multiple deterministic moduli for critical upper HE bit-widths;
- detect tail/off-by-one writes with guard canaries;
- keep the suite practical enough to run after every meaningful code change.

## Fixed regression matrix

### General modulus classes

    20, 30, 40, 50, 52, 55, 59, 60, 61, 62 bits

### Barrett dispatch coverage

`EltwiseMulMod`, `EltwiseFMAMod`, and arbitrary-uint64 `EltwiseReduceMod`
exercise every modulus bit-width:

    20, 21, 22, ..., 61, 62 bits

This maps to all HE-relevant public Barrett dispatch shifts:

    Shift = BitWidth(q) - 2 = 18 ... 60

The SVE white-box arithmetic test additionally instantiates and validates
`ShiftRight128LowPart<Shift>` for **every Shift from 1 through 60**.

For the critical upper widths:

    50, 52, 59, 60, 61, 62 bits

three deterministic NTT-friendly moduli are used: one near the lower quarter of
the bit interval, one near the middle, and one near the upper end. Prime search
is deterministic and does not use random state.

### Element-wise lengths

The matrix includes zero/small lengths, HE-sized vectors through 32768, their
N-1/N/N+1 boundaries, and SVE/unroll boundaries around:

    1*VL, 2*VL, 3*VL, 4*VL, 7*VL, 8*VL, 16*VL

For each of those, VL-1/VL/VL+1 style cases are included where applicable.

### NTT degrees

All selected powers of two from 8 through 32768 are covered:

    8, 16, 32, 64, 128, 256, 512, 1024, 2048,
    4096, 8192, 16384, 32768

The forward factor-4 and inverse factor-2 lazy-input tests run across this full
degree set.

### Deterministic random sampling

The regression uses the fixed seed:

    0x415848454c   // "AXHEL"

There is no seed override in V5. `std::mt19937_64` is combined with a local
rejection-sampling helper rather than `std::uniform_int_distribution`, keeping
the pseudo-random streams stable and independent of distribution implementation.

Main random loops use 256 repetitions; NTT round-trips use 32 repetitions and
the independent O(N^2) polynomial oracle uses 64 repetitions.

## Canary checks

Element-wise Add/Sub/Mul/FMA and Reduce tests place guard words before and after
the payload. The guards are checked after every selected tail size, including
62-bit/factor-4 or factor-8 paths. Input payloads are also checked for accidental
modification when they are not output aliases.

This catches, among other errors:

- SVE tail stores beyond `n`;
- off-by-one writes;
- writes before the output pointer;
- accidental source-buffer modification.

## Public API coverage

- `EltwiseAddMod` vector/vector and vector/scalar
- `EltwiseSubMod` vector/vector and vector/scalar
- `EltwiseMulMod`, factors 1/2/4
- `EltwiseFMAMod`, factors 1/2/4/8 and `op3 == nullptr`
- `EltwiseReduceMod`, identity, 2->1, 4->1, 4->2, arbitrary-uint64 Barrett path
- `NTT::ComputeForward`, normalized/lazy and factor-4 lazy input
- `NTT::ComputeInverse`, normalized/lazy and factor-2 lazy input
- in-place/out-of-place behavior, aliasing, unaligned buffers and n=0

## Internal / white-box coverage

- Native `ReduceInputNative<1/2/4/8>` and reduction helpers
- SVE `ReduceInputSVE<1/2/4/8>` and reduction helpers
- Shoup quotient/lazy multiplication over multiple deterministic q values
- `MulU64ToU128SVE`
- `ShiftRight128LowPart<1>` ... `ShiftRight128LowPart<60>`
- Native forward/inverse Harvey butterflies
- real SVE forward/inverse butterfly implementations through test-only hooks

## Independent integration oracle

The suite retains:

- fixed known NTT vector N=8, q=17, psi=3;
- normalized NTT/INTT round-trip;
- lazy forward -> 4->2 reduction -> inverse round-trip;
- normalized/lazy range invariants;
- NTT -> dyadic multiply -> INTT checked against an independent O(N^2)
  negacyclic convolution modulo `(X^N + 1, q)`.

SVE correctness is therefore not inferred only by comparison with Native.

## Build model

Both test builds are Release (`-O3 -DNDEBUG`) and use the same
`-mcpu=${AXHEL_CPU}` target as the AXHEL library.

- Native regression build: automatic SVE capability detection still runs, then
  `AXHEL_TEST_FORCE_NATIVE=ON` selects the Native backend only for this test build.
- SVE regression build: `AXHEL_TEST_FORCE_NATIVE=OFF`; normal AXHEL automatic SVE
  detection must succeed.

The generated `axhel/util/defines.hpp` is the sole source of truth for
`AXHEL_HAS_SVE`.

## Prerequisite

On Ubuntu:

    sudo apt update
    sudo apt install libgtest-dev

## Configure and build

From the repository root:

    chmod +x tests/scripts/*.sh
    ./tests/scripts/configure-tests.sh

This creates clean builds:

    build-test-native/
    build-test-sve/

On an SVE host the expected configuration is:

Native build:

    AXHEL_SVE_DETECTED: 1
    AXHEL_HAS_SVE (effective backend): 0

SVE build:

    AXHEL_SVE_DETECTED: 1
    AXHEL_HAS_SVE (effective backend): 1

## Run the only test mode

    ./tests/scripts/run-regression.sh

The runner rebuilds both backends before execution and then runs the complete
fixed regression matrix on Native and SVE.

A specific area can still be selected manually, for debugging only:

    ctest --test-dir build-test-sve -R EltwiseMulMod --output-on-failure
    ctest --test-dir build-test-sve -R InternalSVE --output-on-failure
    ctest --test-dir build-test-sve -R NTT --output-on-failure
    ctest --test-dir build-test-sve -R PolynomialMultiply --output-on-failure

## Normal AXHEL builds

With `AXHEL_TESTING=OFF` (default), test-only NTT hooks are not compiled and the
normal AXHEL backend selection remains unchanged.
