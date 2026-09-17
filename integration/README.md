# AXHEL Integrations

This directory contains integration packages for using AXHEL with external homomorphic encryption libraries.

Integrations are organized by target library and version. Each integration is maintained separately from the AXHEL core library and may include version-specific patches and setup instructions.

> [!NOTE]
> Currently, AXHEL provides integration support for Microsoft SEAL versions 4.4.0 through 4.4.5.

## Microsoft SEAL

The Microsoft SEAL integration redirects supported low-level modular arithmetic and Number Theoretic Transform (NTT) operations to AXHEL.

The integration introduces the `SEAL_USE_AXHEL` CMake option and uses AXHEL as an external dependency. AXHEL and Intel HEXL are alternative acceleration backends and cannot be enabled simultaneously.

The integration operates on low-level primitives shared by the homomorphic encryption schemes implemented in Microsoft SEAL and is therefore not specific to a single scheme.

### Supported versions

Each Microsoft SEAL version requires its corresponding integration patch.

| Microsoft SEAL version | AXHEL patch |
|---|---|
| [4.4.0](https://github.com/microsoft/SEAL/releases/tag/v4.4.0) | [`axhel-seal-4.4.0.patch`](seal-4.4.0/axhel-seal-4.4.0.patch) |
| [4.4.1](https://github.com/microsoft/SEAL/releases/tag/v4.4.1) | [`axhel-seal-4.4.1.patch`](seal-4.4.1/axhel-seal-4.4.1.patch) |
| [4.4.2](https://github.com/microsoft/SEAL/releases/tag/v4.4.2) | [`axhel-seal-4.4.2.patch`](seal-4.4.2/axhel-seal-4.4.2.patch) |
| [4.4.3](https://github.com/microsoft/SEAL/releases/tag/v4.4.3) | [`axhel-seal-4.4.3.patch`](seal-4.4.3/axhel-seal-4.4.3.patch) |
| [4.4.4](https://github.com/microsoft/SEAL/releases/tag/v4.4.4) | [`axhel-seal-4.4.4.patch`](seal-4.4.4/axhel-seal-4.4.4.patch) |
| [4.4.5](https://github.com/microsoft/SEAL/releases/tag/v4.4.5) | [`axhel-seal-4.4.5.patch`](seal-4.4.5/axhel-seal-4.4.5.patch) |

> [!IMPORTANT]
> Use the patch matching the exact Microsoft SEAL version being built.

## Installation

The following procedure uses custom installation prefixes for AXHEL and Microsoft SEAL.
For a system-wide installation, the default CMake installation prefix may be used and elevated privileges may be required when running `cmake --install`.

### 1. Build and install AXHEL

Download and extract the AXHEL source archive, then navigate to the AXHEL source directory:

```bash
cd axhel-1.0.0
```

Configure a Release build using the default AXHEL configuration and specify the installation prefix:

```bash
rm -rf build

cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
```

Build and install AXHEL:

```bash
cmake --build build
cmake --install build
```

The default configuration targets the native CPU and automatically enables the SVE backend when supported by the configured target.

### 2. Download Microsoft SEAL

Download and extract the source archive of one of the supported Microsoft SEAL versions.

Navigate to the Microsoft SEAL source directory. For example, for Microsoft SEAL 4.4.0:

```bash
cd SEAL-4.4.0
```

### 3. Verify the integration patch

Before applying the patch, verify that it is compatible with the selected Microsoft SEAL source tree:

```bash
git apply --check /path/to/axhel-1.0.0/integration/seal-4.4.0/axhel-seal-4.4.0.patch
```

If the command completes without errors, the patch can be applied.

### 4. Apply the integration patch

Apply the version-specific patch:

```bash
git apply /path/to/axhel-1.0.0/integration/seal-4.4.0/axhel-seal-4.4.0.patch
```
> [!NOTE]
> For a different supported Microsoft SEAL release, replace both the source directory and patch path with the corresponding version.

### 5. Configure Microsoft SEAL with AXHEL

Configure Microsoft SEAL with AXHEL enabled:

```bash
rm -rf build-axhel

cmake -S . -B build-axhel -DSEAL_USE_AXHEL=ON -DSEAL_USE_INTEL_HEXL=OFF -DCMAKE_PREFIX_PATH=/axhel/installation/path -DCMAKE_INSTALL_PREFIX=/seal/path/to/install
```

For AXHEL installations in non-standard locations, `CMAKE_PREFIX_PATH` must include the AXHEL installation prefix so that Microsoft SEAL can locate the AXHEL CMake package configuration files.

`SEAL_USE_AXHEL=ON` enables the AXHEL integration, while `SEAL_USE_INTEL_HEXL=OFF` ensures that the two alternative acceleration backends are not enabled simultaneously.

### 6. Build and install Microsoft SEAL

Build the patched Microsoft SEAL library:

```bash
cmake --build build-axhel
```

Install Microsoft SEAL:

```bash
cmake --install build-axhel
```

The resulting Microsoft SEAL installation is configured to use the AXHEL backend for the supported modular arithmetic and NTT operations.

## Using the installed Microsoft SEAL package

When the AXHEL-enabled Microsoft SEAL installation is used by another CMake project, both the Microsoft SEAL and AXHEL installation prefixes must be discoverable by CMake.

For example:

```bash
cmake -S . -B build -DSEAL_ROOT=/seal/installation/path -DCMAKE_PREFIX_PATH=/axhel/installation/path
```

This allows the Microsoft SEAL package configuration to locate its AXHEL dependency.