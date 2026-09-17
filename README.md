![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![ARM](https://img.shields.io/badge/Architecture-AArch64-green.svg)
![SVE](https://img.shields.io/badge/ARM-SVE-success.svg)

# ARM aXceleration for Homomorphic Encryption Library

AXHEL is an open-source C++ acceleration library for AArch64 processors providing optimized modular arithmetic and Number Theoretic Transform (NTT) kernels for homomorphic encryption.

AXHEL exploits the Arm Scalable Vector Extension (SVE) when available while providing portable native/scalar implementations as fallback. It is designed as a lightweight acceleration layer that can be integrated into existing homomorphic encryption frameworks.

## Introduction

Homomorphic Encryption (HE) enables computations to be performed directly on encrypted data without revealing the underlying cleartext. 

Modern lattice-based HE schemes, such as CKKS, BFV and BGV, operate on polynomials over the quotient ring $\mathbb{Z}_q[X]/(X^N+1)$. 
In practical implementations, the coefficient modulus $q$ is represented as a chain of smaller prime moduli using the Residue Number System (RNS). 
Each polynomial coefficient is therefore represented as a set of independent residues, typically fitting within native 64-bit unsigned integers, 
with individual moduli usually ranging from 30 to 60 bits.

Although ciphertext operations are expressed as polynomial additions, multiplications and transformations, they ultimately reduce to a large number 
of modular arithmetic operations on vectors of 64-bit coefficients. Since these operations are executed repeatedly throughout the evaluation process, 
they represent one of the primary computational bottlenecks of practical HE implementations.

AXHEL (Arm aXceleration for Homomorphic Encryption Library) is an open-source C++ library that provides optimized implementations of low-level modular arithmetic and Number Theoretic Transform (NTT) kernels for Arm AArch64 processors. By exploiting Arm Scalable Vector Extension (SVE), AXHEL accelerates these computational building blocks of higher-level homomorphic encryption operations, while preserving portable scalar implementations for maximum compatibility.

These kernels are not specific to a single homomorphic encryption scheme. They are shared by schemes such as CKKS, BFV, and BGV, allowing AXHEL to accelerate their execution whenever the corresponding low-level arithmetic and NTT primitives are used.

## Features

The library provides the following kernels:

| Kernel | Native | SVE | Description |
|---|:---:|:---:|---|
| `EltwiseAddMod` | ✓ | ✓ | Element-wise modular addition |
| `EltwiseSubMod` | ✓ | ✓ | Element-wise modular subtraction |
| `EltwiseMulMod` | ✓ | ✓ | Element-wise modular multiplication |
| `EltwiseFMAMod` | ✓ | ✓ | Element-wise modular multiply-add |
| `EltwiseReduceMod` | ✓ | ✓ | Element-wise modular reduction |
| `ForwardNTT` | ✓ | ✓ | Forward negacyclic Number Theoretic Transform |
| `InverseNTT` | ✓ | ✓ | Inverse negacyclic Number Theoretic Transform |

Addition and subtraction support vector-vector and vector-scalar forms; multiplication is vector-vector, while FMA performs vector-scalar multiplication with an optional vector addend. The NTT interface provides lazy and normalized forward and inverse negacyclic transforms.

Each kernel is available through native and, when supported, SVE-optimized implementations sharing the same public API. During CMake configuration, AXHEL detects SVE support for the configured target and enables the SVE backend when available; otherwise, the native implementation is used automatically.

## Building AXHEL

AXHEL uses CMake as its build system, enabling a portable and configurable build process across supported Arm AArch64 platforms. 
The following sections describe the required dependencies, available build options, and the steps needed to configure, compile and install the library.

### Requirements

The following software is required to build AXHEL.

| Dependency | Requirement |
|------------|-------------|
| CMake | 3.13 or later |
| Compiler | GCC 10+ or Clang 12+ |

AXHEL has been developed and tested on Linux-based Arm AArch64 platforms. Support for Arm Scalable Vector Extension (SVE) is automatically detected during CMake configuration. 
When SVE is supported by the configured target, the optimized SVE backend is enabled; otherwise, AXHEL falls back to the portable native implementation.

### Compile-time options

AXHEL supports the following library-specific compile-time options.

| AXHEL Option | Values | Default | Description |
|--------------|--------|:-------:|-------------|
| `AXHEL_SHARED_LIB` | `ON`, `OFF` | `OFF` | Build AXHEL as a shared library instead of a static library. |
| `AXHEL_TREAT_WARNING_AS_ERROR` | `ON`, `OFF` | `OFF` | Treat compiler warnings as errors. |
| `AXHEL_CPU` | Any CPU name supported by the compiler, e.g. `native`, `neoverse-v1`, `neoverse-v2` | `native` | Select the target ARM CPU passed to the compiler through `-mcpu=<value>`. |
| `AXHEL_OPT_REPORT` | `ON`, `OFF` | `OFF` | Enable compiler optimization reports. |

In addition, AXHEL supports the standard CMake configuration variables.

| CMake Variable | Typical Values | Default | Description |
|----------------|----------------|:-------:|-------------|
| `CMAKE_BUILD_TYPE` | `Release`, `Debug` | `Release` | Select the build configuration. `Debug` enables debug symbols, AXHEL runtime tracing, and AddressSanitizer. |
| `CMAKE_INSTALL_PREFIX` | `</install/path>` | System default | Specify the installation directory used by `cmake --install`. |


### Compile AXHEL
AXHEL supports building from source using the CMake build system. The following instructions describe how to configure, build and install the library.

1. **Navigate to the project root directory.**

    After cloning or downloading the repository:

    ```bash
    cd axhel
    ```

2. **Configure the library.**

    ```bash
    cmake -S . -B build
    ```
    This creates a Release build targeting the native CPU with the default AXHEL configuration.

    Additional compile-time options can be specified using `-D`. 

    For example, to install AXHEL in a custom location, configure the build with:

    ```bash
    cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
    ```

    For example, to build AXHEL for a Neoverse V1 Arm CPU:

    ```bash
    cmake -S . -B build -DAXHEL_CPU=neoverse-v1
    ```

3. **Build the library.**

    ```bash
    cmake --build build
    ```

    This command builds the AXHEL library and any enabled targets in the `build/` directory.

4. **Install the library.**

    ```bash
    cmake --install build
    ```

    The installation includes the AXHEL library, public headers, CMake package configuration files, and exported CMake targets, 
    allowing AXHEL to be easily integrated into external CMake projects using `find_package(AXHEL)`.

## Example

A standalone example demonstrating the main AXHEL operations is provided in the [`example/`](example/) directory.

The example covers element-wise modular arithmetic, modular reduction, and forward/inverse NTT operations, including lazy and normalized transforms.

Navigate to the example directory:

```bash
cd axhel/example
```

Configure and build the example:

```bash
cmake -S . -B build
cmake --build build
```

Run the example:

```bash
./build/example
```

AXHEL is built automatically as part of the example build and does not need to be installed separately.

## Testing

AXHEL includes a deterministic regression suite covering both the Native and SVE backends. The tests validate the public modular arithmetic and NTT APIs, internal arithmetic helpers, boundary conditions, aliasing, and SVE-specific implementations.

The regression suite requires [GoogleTest](https://github.com/google/googletest). 

From the AXHEL root folder, run the complete regression suite with:

```bash
chmod +x tests/scripts/*.sh
./tests/scripts/run-regression.sh
```

The script builds and tests both the Native and SVE backends using separate Release configurations.

> [!NOTE]
> The SVE regression requires an SVE-capable target. For the complete regression matrix, configuration options, and advanced usage, see the [tests README](tests/README.md).


## Debugging

For maximum performance, AXHEL performs only minimal runtime validation in `Release` builds. 
To debug AXHEL, configure and build the library with `-DCMAKE_BUILD_TYPE=Debug`. 

```bash
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug
```

This generates a debug version of the library (e.g., `libaxhel_debug.a`) with debug symbols, enables internal `AXHEL_LOG` tracing, and links against AddressSanitizer.

Enabling `CMAKE_BUILD_TYPE=Debug` introduces a significant runtime overhead and is intended exclusively for debugging and development, not for performance evaluation.

## Integration

AXHEL is designed to be integrated into homomorphic encryption frameworks requiring high-performance modular arithmetic and NTT kernels.

The `integration/` directory provides version-specific integration files for supported external libraries. Each integration is kept separate from the AXHEL
core library and can be applied after AXHEL has been built and installed.

The complete list of supported integrations, versions, and setup instructions is maintained in the [integration README](integration/README.md).

> [!NOTE]
> Currently, AXHEL provides integration support for Microsoft SEAL versions 4.4.0 through 4.4.5.



## License

AXHEL is distributed under the Apache License 2.0.

Some source files are derived from or inspired by Intel HEXL and retain the corresponding copyright notices in accordance with the Apache License.

All Arm-specific implementations and additional developments are Copyright © 2026 University of Pisa.

