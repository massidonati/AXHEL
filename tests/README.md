# AXHEL Tests

AXHEL includes a deterministic regression suite for validating the correctness of the Native and SVE implementations.

The tests cover the public modular arithmetic and Number Theoretic Transform (NTT) APIs together with selected internal arithmetic primitives used by the optimized backends. The suite is intended to detect functional regressions after changes to the AXHEL implementation.

## Test coverage

The regression suite covers:

- element-wise modular addition and subtraction, including vector-vector and vector-scalar forms;
- modular multiplication and fused multiply-add;
- modular reduction for the supported input and output ranges;
- forward and inverse NTT operations, including lazy and normalized transforms;
- in-place and out-of-place execution;
- aliasing, unaligned buffers, boundary lengths, and zero-length inputs;
- Native arithmetic and NTT helpers;
- SVE arithmetic helpers and NTT implementations when SVE is available;
- range and round-trip properties of the NTT;
- polynomial multiplication checked against an independent negacyclic convolution reference.

The tests use deterministic input generation so that repeated executions exercise the same test cases.

Guard values are also used around selected buffers to detect out-of-bounds and off-by-one writes.

## Requirements

The test suite uses [GoogleTest](https://github.com/google/googletest).

For example, On Ubuntu it can be installed with:

```bash
sudo apt update
sudo apt install libgtest-dev
```

Alternatively, GoogleTest can be downloaded automatically by CMake by configuring AXHEL with:

```bash
-DAXHEL_TEST_FETCH_GTEST=ON
```

## Test configuration options

The following CMake options control the AXHEL test build:

| Option | Default | Description |
|---|:---:|---|
| `AXHEL_TESTING` | `OFF` | Enables the AXHEL correctness and regression tests |
| `AXHEL_TEST_FETCH_GTEST` | `OFF` | Downloads GoogleTest if it is not available on the system |
| `AXHEL_TEST_FORCE_NATIVE` | `OFF` | Forces the Native backend for a test build |
| `AXHEL_CPU` | `native` | Selects the Arm CPU target passed to the compiler through `-mcpu` |

`AXHEL_TEST_FORCE_NATIVE` is valid only when `AXHEL_TESTING=ON`.

Normal AXHEL builds are unaffected because testing is disabled by default.

## Test scripts

Two helper scripts are provided in [`tests/scripts/`](scripts/):

- `configure-tests.sh` configures and builds separate Native and SVE test configurations;
- `run-regression.sh` builds both configurations and executes the complete regression suite.

### Script defaults

The scripts use the following defaults:

| Setting | Default |
|---|---|
| C++ compiler | `g++` |
| CPU target | `native` |
| Build type | `Release` |
| Native build directory | `build-test-native/` |
| SVE build directory | `build-test-sve/` |
| Parallel build jobs | CMake default parallelism |

The compiler, CPU target, and number of parallel build jobs can be changed through environment variables:

| Environment variable | Description |
|---|---|
| `CXX` | C++ compiler used to build AXHEL and the tests |
| `AXHEL_CPU` | CPU target passed to AXHEL through `-mcpu` |
| `AXHEL_TEST_JOBS` | Number of parallel build jobs |

For example:

```bash
CXX=g++ AXHEL_CPU=neoverse-v2 AXHEL_TEST_JOBS=8 \
./tests/scripts/configure-tests.sh
```

## Configure and build

From the AXHEL repository root:

```bash
chmod +x tests/scripts/*.sh
./tests/scripts/configure-tests.sh
```

The script creates two clean build directories:

```text
build-test-native/
build-test-sve/
```

Both configurations use:

```text
CMAKE_BUILD_TYPE=Release
AXHEL_TESTING=ON
```

For the Native build, the script sets:

```text
AXHEL_TEST_FORCE_NATIVE=ON
```

while the SVE build uses the normal AXHEL SVE detection with:

```text
AXHEL_TEST_FORCE_NATIVE=OFF
```

The script verifies that the effective backend matches the requested test configuration. The SVE build therefore requires an SVE-capable compiler and target.

## Run the regression suite

From the repository root:

```bash
./tests/scripts/run-regression.sh
```

If the Native or SVE build directory does not exist, the runner automatically invokes `configure-tests.sh` first.

Before running the tests, both configurations are rebuilt to ensure that the regression suite does not use stale binaries.

The complete regression suite is then executed for both the Native and SVE backends.

## Run individual tests

Individual test areas can be selected with CTest when debugging a specific part of AXHEL.

Use:

```bash
ctest --test-dir build-test-sve -R <area> --output-on-failure
```

where `<area>` can be one of: `EltwiseAddMod`, `EltwiseSubMod`, `EltwiseAddSubMod`, `EltwiseMulMod`, `EltwiseFMAMod`, `EltwiseOperations`, `EltwiseReduceMod`, `NTT`, `PolynomialMultiply`, `InternalNativeReduction`, `InternalShoupNative`, `InternalNTTButterflyNative`, `InternalSVEArithmetic`, `InternalSVEReduction`, `InternalSVEShoup`, or `InternalSVEButterfly`.

For example:

```bash
ctest --test-dir build-test-sve -R EltwiseMulMod --output-on-failure
```

For the Native build, replace `build-test-sve` with `build-test-native`.

## Native-only testing

On a system where SVE is not available, a Native-only test build can be configured directly with CMake:

```bash
cmake -S . -B build-test-native \
    -DCMAKE_BUILD_TYPE=Release \
    -DAXHEL_TESTING=ON \
    -DAXHEL_TEST_FORCE_NATIVE=ON
```

Build and run the tests with:

```bash
cmake --build build-test-native --parallel
ctest --test-dir build-test-native -L regression --output-on-failure
```

## Normal AXHEL builds

Testing is disabled by default:

```text
AXHEL_TESTING=OFF
```

When testing is disabled, test-only sources and NTT hooks are not compiled and the normal AXHEL backend selection remains unchanged.
