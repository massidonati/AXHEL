#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
JOBS="${AXHEL_TEST_JOBS:-}"
NATIVE_BUILD="$ROOT/build-test-native"
SVE_BUILD="$ROOT/build-test-sve"

if [[ ! -d "$NATIVE_BUILD" || ! -d "$SVE_BUILD" ]]; then
    "$ROOT/tests/scripts/configure-tests.sh"
fi

build_one() {
    local build="$1"
    if [[ -n "$JOBS" ]]; then
        cmake --build "$build" --parallel "$JOBS"
    else
        cmake --build "$build" --parallel
    fi
}

# Always rebuild first: the regression gate must never exercise stale binaries.
build_one "$NATIVE_BUILD"
build_one "$SVE_BUILD"

echo "=== Native extended regression ==="
ctest --test-dir "$NATIVE_BUILD" -L regression --output-on-failure

echo
echo "=== SVE extended regression (public + internal) ==="
ctest --test-dir "$SVE_BUILD" -L regression --output-on-failure

echo
echo "PASS: deterministic extended Native and SVE regression."
