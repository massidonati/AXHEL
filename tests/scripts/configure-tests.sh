#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CXX="${CXX:-g++}"
CPU="${AXHEL_CPU:-native}"
JOBS="${AXHEL_TEST_JOBS:-}"

build_cmd() {
    local build="$1"
    if [[ -n "$JOBS" ]]; then
        cmake --build "$build" --parallel "$JOBS"
    else
        cmake --build "$build" --parallel
    fi
}

verify_backend() {
    local impl="$1"
    local build="$2"
    local defines="$build/axhel/include/axhel/util/defines.hpp"

    if [[ ! -f "$defines" ]]; then
        echo "ERROR: generated AXHEL configuration header not found: $defines" >&2
        exit 1
    fi

    if [[ "$impl" == "sve" ]]; then
        if ! grep -q '^#define AXHEL_HAS_SVE' "$defines"; then
            echo "ERROR: SVE test build requested, but AXHEL automatic detection did not enable SVE." >&2
            echo "       Check compiler/CPU support and AXHEL_CPU=${CPU}." >&2
            exit 1
        fi
    else
        if grep -q '^#define AXHEL_HAS_SVE' "$defines"; then
            echo "ERROR: forced-Native test build unexpectedly has AXHEL_HAS_SVE enabled." >&2
            exit 1
        fi
    fi
}

configure() {
    local impl="$1"
    local build="$ROOT/build-test-$impl"

    # A clean configure avoids stale generated defines or source lists after
    # changing backend/test configuration.
    rm -rf "$build"

    local force_native=OFF
    if [[ "$impl" == "native" ]]; then
        force_native=ON
    fi

    echo
    echo "=== Configure ${impl} test build ==="
    cmake -S "$ROOT" -B "$build" \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$CXX" \
      -DAXHEL_CPU="$CPU" \
      -DAXHEL_TESTING=ON \
      -DAXHEL_TEST_FORCE_NATIVE="$force_native"

    verify_backend "$impl" "$build"
    build_cmd "$build"
}

configure native
configure sve

echo
echo "PASS: configured and built build-test-native and build-test-sve"
