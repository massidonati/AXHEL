#!/usr/bin/env bash
# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REPETITIONS=100
WARMUP=10
BATCH=1
PMU_PROFILE="none"
PMU_CNTVCT=0
CPU="auto"
DEGREES="2048,4096,8192,16384,32768"
BITS="20,25,30,40,50,52,54,56,58,60,61,62"
KERNELS="all"
AXHEL_CPU_TARGET="native"
CXX_COMPILER="${CXX:-g++}"
JOBS="${AXHEL_BENCH_JOBS:-}"
STAMP="$(date +%Y%m%d-%H%M%S)"
OUTPUT_DIR="${ROOT_DIR}/benchmark/results/${STAMP}"
WRITE_RAW=1
SHUFFLE_SCENARIOS=0
SHUFFLE_SEED=0

usage() {
    cat <<USAGE
Usage: $0 [options]
  --repetitions N
  --warmup N
  --batch N                    kernel invocations per PMU sample (default 1)
  --pmu-profile P              none|l1|l2|stalls (default none)
  --cntvct                     read CNTVCT_EL0 in diagnostic profiles
  --cpu auto|off|ID
  --degrees LIST
  --bits LIST
  --kernels LIST|all
  --axhel-cpu CPU_TARGET       value passed to -DAXHEL_CPU= (default native)
  --shuffle-scenarios         randomize (degree, qi_bits) scenario order
  --shuffle-seed N            deterministic scenario-order seed
  --output DIR
  --no-raw
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repetitions) REPETITIONS="$2"; shift 2 ;;
        --warmup) WARMUP="$2"; shift 2 ;;
        --batch) BATCH="$2"; shift 2 ;;
        --pmu-profile) PMU_PROFILE="$2"; shift 2 ;;
        --cntvct) PMU_CNTVCT=1; shift ;;
        --cpu) CPU="$2"; shift 2 ;;
        --degrees) DEGREES="$2"; shift 2 ;;
        --bits) BITS="$2"; shift 2 ;;
        --kernels) KERNELS="$2"; shift 2 ;;
        --axhel-cpu) AXHEL_CPU_TARGET="$2"; shift 2 ;;
        --shuffle-scenarios) SHUFFLE_SCENARIOS=1; shift ;;
        --shuffle-seed) SHUFFLE_SEED="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --no-raw) WRITE_RAW=0; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

mkdir -p "$OUTPUT_DIR/scalar-strict" "$OUTPUT_DIR/native" "$OUTPUT_DIR/sve"

build_target() {
    local build_dir="$1"
    if [[ -n "$JOBS" ]]; then
        cmake --build "$build_dir" --target axhel-benchmark --parallel "$JOBS"
    else
        cmake --build "$build_dir" --target axhel-benchmark --parallel
    fi
}

verify_backend() {
    local expected="$1"
    local build_dir="$2"
    local defines="$build_dir/axhel/include/axhel/util/defines.hpp"
    if [[ ! -f "$defines" ]]; then
        echo "ERROR: generated AXHEL configuration header not found: $defines" >&2
        exit 1
    fi
    if [[ "$expected" == "sve" ]]; then
        if ! grep -q '^#define AXHEL_HAS_SVE' "$defines"; then
            echo "ERROR: SVE benchmark build requested, but AXHEL automatic detection did not enable SVE." >&2
            echo "       Check compiler/CPU support and --axhel-cpu $AXHEL_CPU_TARGET." >&2
            exit 1
        fi
    else
        if grep -q '^#define AXHEL_HAS_SVE' "$defines"; then
            echo "ERROR: $expected benchmark build unexpectedly has AXHEL_HAS_SVE enabled." >&2
            exit 1
        fi
    fi
}

verify_scalar_compile_flags() {
    local build_dir="$1"
    local cc="$build_dir/compile_commands.json"
    if [[ ! -f "$cc" ]]; then
        echo "ERROR: compile_commands.json not found for scalar verification: $cc" >&2
        exit 1
    fi

    # Verify the flags are attached to AXHEL library translation units, not
    # merely to the benchmark driver executable.
    local axhel_commands
    axhel_commands="$(grep -E 'axhel/(eltwise|ntt|number-theory)/.*\.cpp' "$cc" || true)"
    if [[ -z "$axhel_commands" ]]; then
        echo "ERROR: no AXHEL library compile commands found for scalar verification." >&2
        exit 1
    fi
    if echo "$axhel_commands" | grep -E '(-sve\.cpp|ntt-sve\.cpp)' >/dev/null; then
        echo "ERROR: strict-scalar build unexpectedly contains explicit SVE translation units." >&2
        exit 1
    fi

    if "$CXX_COMPILER" --version 2>/dev/null | head -1 | grep -qi clang; then
        if ! echo "$axhel_commands" | grep -- '-fno-vectorize' >/dev/null || \
           ! echo "$axhel_commands" | grep -- '-fno-slp-vectorize' >/dev/null; then
            echo "ERROR: strict-scalar Clang vectorization-disable flags not found on AXHEL sources." >&2
            exit 1
        fi
    else
        if ! echo "$axhel_commands" | grep -- '-fno-tree-loop-vectorize' >/dev/null || \
           ! echo "$axhel_commands" | grep -- '-fno-tree-slp-vectorize' >/dev/null; then
            echo "ERROR: strict-scalar GCC vectorization-disable flags not found on AXHEL sources." >&2
            exit 1
        fi
    fi
}

COMMON_CMAKE=(
    -DCMAKE_BUILD_TYPE=Release
    -DAXHEL_BENCHMARK=ON
    -DAXHEL_TESTING=OFF
    -DAXHEL_CPU="$AXHEL_CPU_TARGET"
    -DCMAKE_CXX_COMPILER="$CXX_COMPILER"
)

COMMON_RUN=(
    --repetitions "$REPETITIONS"
    --warmup "$WARMUP"
    --batch "$BATCH"
    --pmu-profile "$PMU_PROFILE"
    --cpu "$CPU"
    --degrees "$DEGREES"
    --bits "$BITS"
    --kernels "$KERNELS"
)
if [[ "$SHUFFLE_SCENARIOS" -eq 1 ]]; then
    COMMON_RUN+=(--shuffle-scenarios --shuffle-seed "$SHUFFLE_SEED")
fi
if [[ "$WRITE_RAW" -eq 0 ]]; then
    COMMON_RUN+=(--no-raw)
fi
if [[ "$PMU_CNTVCT" -eq 1 ]]; then
    COMMON_RUN+=(--cntvct)
fi

echo "== Configure/build ScalarStrict public backend (compiler auto-vectorization OFF) =="
rm -rf "$ROOT_DIR/build-benchmark-scalar-strict"
cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build-benchmark-scalar-strict" \
    "${COMMON_CMAKE[@]}" \
    -DAXHEL_BENCHMARK_FORCE_NATIVE=OFF \
    -DAXHEL_BENCHMARK_FORCE_SCALAR=ON
verify_backend scalar-strict "$ROOT_DIR/build-benchmark-scalar-strict"
verify_scalar_compile_flags "$ROOT_DIR/build-benchmark-scalar-strict"
build_target "$ROOT_DIR/build-benchmark-scalar-strict"

echo "== Run ScalarStrict =="
"$ROOT_DIR/build-benchmark-scalar-strict/benchmark/axhel-benchmark" \
    "${COMMON_RUN[@]}" \
    --expect-backend ScalarStrict \
    --output-dir "$OUTPUT_DIR/scalar-strict"

echo "== Configure/build Native public backend (normal compiler optimization policy) =="
rm -rf "$ROOT_DIR/build-benchmark-native"
cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build-benchmark-native" \
    "${COMMON_CMAKE[@]}" \
    -DAXHEL_BENCHMARK_FORCE_NATIVE=ON \
    -DAXHEL_BENCHMARK_FORCE_SCALAR=OFF
verify_backend native "$ROOT_DIR/build-benchmark-native"
build_target "$ROOT_DIR/build-benchmark-native"

echo "== Run Native =="
"$ROOT_DIR/build-benchmark-native/benchmark/axhel-benchmark" \
    "${COMMON_RUN[@]}" \
    --expect-backend Native \
    --output-dir "$OUTPUT_DIR/native"

echo "== Configure/build SVE public backend =="
rm -rf "$ROOT_DIR/build-benchmark-sve"
cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build-benchmark-sve" \
    "${COMMON_CMAKE[@]}" \
    -DAXHEL_BENCHMARK_FORCE_NATIVE=OFF \
    -DAXHEL_BENCHMARK_FORCE_SCALAR=OFF
verify_backend sve "$ROOT_DIR/build-benchmark-sve"
build_target "$ROOT_DIR/build-benchmark-sve"

echo "== Run SVE =="
"$ROOT_DIR/build-benchmark-sve/benchmark/axhel-benchmark" \
    "${COMMON_RUN[@]}" \
    --expect-backend SVE \
    --output-dir "$OUTPUT_DIR/sve"

echo "== Compute ScalarStrict/Native/SVE speedups =="
python3 "$ROOT_DIR/benchmark/scripts/compare_results.py" \
    "$OUTPUT_DIR/scalar-strict/summary.csv" \
    "$OUTPUT_DIR/native/summary.csv" \
    "$OUTPUT_DIR/sve/summary.csv" \
    "$OUTPUT_DIR/speedup.csv"

cat > "$OUTPUT_DIR/run-config.txt" <<CFG
repetitions=$REPETITIONS
warmup=$WARMUP
batch=$BATCH
pmu_profile=$PMU_PROFILE
cntvct=$PMU_CNTVCT
cpu=$CPU
degrees=$DEGREES
bits=$BITS
kernels=$KERNELS
axhel_cpu=$AXHEL_CPU_TARGET
scenario_order=$([[ "$SHUFFLE_SCENARIOS" -eq 1 ]] && echo randomized || echo fixed)
shuffle_seed=$SHUFFLE_SEED
raw_samples=$WRITE_RAW
baselines=ScalarStrict,Native,SVE
primary_speedup=ScalarStrict/SVE
secondary_speedups=ScalarStrict/Native,Native/SVE
CFG

echo
echo "Results: $OUTPUT_DIR"
echo "  ScalarStrict summary: $OUTPUT_DIR/scalar-strict/summary.csv"
echo "  Native summary:       $OUTPUT_DIR/native/summary.csv"
echo "  SVE summary:          $OUTPUT_DIR/sve/summary.csv"
echo "  Combined speedups:    $OUTPUT_DIR/speedup.csv"
