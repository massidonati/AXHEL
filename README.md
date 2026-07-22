![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)
![ARM](https://img.shields.io/badge/Architecture-AArch64-green.svg)
![SVE](https://img.shields.io/badge/ARM-SVE-success.svg)

# ARM aXceleration for Homomorphic Encryption Library

AXHEL is an open-source C++ library providing optimized modular arithmetic kernels for homomorphic encryption on ARM processors. AXHEL accelerates the arithmetic primitives commonly used by homomorphic encryption libraries by exploiting ARM Scalable Vector Extension (SVE) while maintaining portable scalar implementations.

AXHEL is designed as a lightweight acceleration layer that can be integrated into existing homomorphic encryption frameworks to improve the performance of modular arithmetic on modern AArch64 platforms.

## Introduction

Homomorphic Encryption (HE) enables computations to be performed directly on encrypted data without revealing the underlying cleartext. 

Modern lattice-based HE schemes, such as CKKS, BFV and BGV, operate on polynomials over the quotient ring $\mathbb{Z}_q[X]/(X^N+1)$. 
In practical implementations, the coefficient modulus $q$ is represented as a chain of smaller prime moduli using the Residue Number System (RNS). 
Each polynomial coefficient is therefore represented as a set of independent residues, typically fitting within native 64-bit unsigned integers, 
with individual moduli usually ranging from 30 to 60 bits.

Although ciphertext operations are expressed as polynomial additions, multiplications and transformations, they ultimately reduce to a large number 
of modular arithmetic operations on vectors of 64-bit coefficients. Since these operations are executed repeatedly throughout the evaluation process, 
they represent one of the primary computational bottlenecks of practical HE implementations.

AXHEL (ARM aXceleration for Homomorphic Encryption Library) is an open-source C++ library that provides optimized implementations of these arithmetic kernels 
for ARM AArch64 processors. By exploiting ARM Scalable Vector Extension (SVE), AXHEL accelerates the low-level modular arithmetic primitives that constitute 
the building blocks of higher-level homomorphic encryption operations, while preserving portable scalar implementations for maximum compatibility.

## Features

The library covers the following kernels:
| Kernel | Scalar | SVE |
|---|---:|---:|
| EltwiseAddMod | ✓ | ✓ |
| EltwiseSubMod | ✓ | ✓ |
| EltwiseMulMod | ✓ | ✓ |
| EltwiseFMAMod | ✓ | ✓ |
| EltwiseReduceMod | ✓ | ✓ |
| NTT | ✓ | ✓ |
| INTT | ✓ | ✓ |

Each kernel is available through multiple implementations sharing the same public API. During the CMake configuration phase, AXHEL automatically detects the capabilities 
of the target compiler and processor and selects the most appropriate implementation. This approach enables architecture-specific optimizations while preserving portability. 

## Building AXHEL

AXHEL uses CMake as its build system, enabling a portable and configurable build process across supported ARM AArch64 platforms. 
The following sections describe the required dependencies, available build options, and the steps needed to configure, compile and install the library.

### Requirements

The following software is required to build AXHEL.

| Dependency | Requirement |
|------------|-------------|
| CMake | 3.13 or later |
| Compiler | GCC 12+ or Clang 15+ |

AXHEL has been developed and tested on Linux-based ARM AArch64 platforms. Support for ARM Scalable Vector Extension (SVE) is automatically detected during the CMake configuration process. 
When supported by the target compiler and processor, the corresponding optimized kernels are enabled automatically; otherwise, AXHEL transparently falls back to the portable scalar implementation.

### Compile-time options

AXHEL supports the following library-specific compile-time options.

| AXHEL Option | Values | Default | Description |
|--------------|--------|:-------:|-------------|
| `AXHEL_SHARED_LIB` | `ON`, `OFF` | `OFF` | Build AXHEL as a shared library instead of a static library. |
| `AXHEL_TREAT_WARNING_AS_ERROR` | `ON`, `OFF` | `OFF` | Treat compiler warnings as errors. |
| `AXHEL_CPU` | Any CPU name supported by the compiler, e.g. `native`, `neoverse-v1`, `neoverse-v2` | `native` | Select the target ARM CPU passed to the compiler through `-mcpu=<value>`. |

In addition, AXHEL supports the standard CMake configuration variables.

| CMake Variable | Typical Values | Default | Description |
|----------------|----------------|:-------:|-------------|
| `CMAKE_BUILD_TYPE` | `Release`, `Debug` | `Release` | Select the build configuration. |
| `CMAKE_INSTALL_PREFIX` | `<path>` | System default | Specify the installation directory used by `cmake --install`. |


### Compile AXHEL
AXHEL supports building from source using the CMake build system. The following instructions describe how to configure, build and install the library.

After cloning or downloading the repository, navigate to the project root directory.

```bash
cd AXHEL
```
Configure the project using CMake. Additional compile-time options can be specified during the CMake configuration step by adding the desired option with the `-D` flag.

```bash
cmake -S . -B build
```
For example, to use a non-default installation directory, configure the build with:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
```

For example, to build AXHEL for the a Neoverse V1 ARM CPU:

```bash
cmake -S . -B build -DAXHEL_CPU=neoverse-v1
```

Build AXHEL by running:

```bash
cmake --build build
```

This command builds the AXHEL library and any enabled targets in the `build/` directory.

Install AXHEL by running:

```bash
cmake --install build
```

The installation includes the AXHEL library, public headers, CMake package configuration files, and exported CMake targets, 
enabling AXHEL to be easily integrated into external CMake projects using `find_package(AXHEL)`.


# Integration

AXHEL is designed to be easily integrated into homomorphic encryption frameworks requiring high-performance modular arithmetic kernels.


# Contributing

Contributions, bug reports and feature requests are welcome.


# License

AXHEL is distributed under the Apache License 2.0.

Some source files are derived from or inspired by Intel HEXL and retain the corresponding copyright notices in accordance with the Apache License.

All ARM-specific implementations and additional developments are Copyright © 2026 University of Pisa.
