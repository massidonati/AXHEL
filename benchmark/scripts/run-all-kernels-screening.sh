#!/usr/bin/env bash

# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

RUNS="${1:-10}"
BASE_SHUFFLE_SEED="${2:-2026082900}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT_DIR}"

TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
CAMPAIGN_DIR="benchmark/results/${TIMESTAMP}-all-kernels-screening-randomized"

mkdir -p "${CAMPAIGN_DIR}"

echo "AXHEL all-kernels randomized screening campaign"
echo "Runs              : ${RUNS}"
echo "Base shuffle seed : ${BASE_SHUFFLE_SEED}"
echo "CPU               : 17"
echo "Repetitions       : 200"
echo "Warmup            : 50"
echo "Batch             : 1"
echo "PMU profile       : none"
echo "Degrees           : 2048,4096,8192,16384,32768"
echo "qi bits           : 20,25,30,40,50,52,54,56,58,60,61,62"
echo "Kernels           : all"
echo "AXHEL CPU target   : native"
echo "Scenario order    : randomized"
echo "Output            : ${CAMPAIGN_DIR}"
echo "Git commit        : $(git rev-parse HEAD)"
echo

# Snapshot the exact source state used for the campaign.
git rev-parse HEAD > "${CAMPAIGN_DIR}/git-commit.txt"
git status --short > "${CAMPAIGN_DIR}/git-status.txt"
git diff > "${CAMPAIGN_DIR}/git-diff.patch"

cat > "${CAMPAIGN_DIR}/campaign-config.txt" <<EOF
runs=${RUNS}
base_shuffle_seed=${BASE_SHUFFLE_SEED}
cpu=17
repetitions=200
warmup=50
batch=1
pmu_profile=none
degrees=2048,4096,8192,16384,32768
bits=20,25,30,40,50,52,54,56,58,60,61,62
kernels=all
axhel_cpu=native
scenario_order=randomized
raw_samples=enabled
EOF

printf "run,shuffle_seed\n" > "${CAMPAIGN_DIR}/shuffle-seeds.csv"

for RUN in $(seq 1 "${RUNS}")
do
    RUN_ID="$(printf "%02d" "${RUN}")"
    SHUFFLE_SEED=$((BASE_SHUFFLE_SEED + RUN))
    OUT="${CAMPAIGN_DIR}/run-${RUN_ID}"

    printf "%s,%s\n" \
        "${RUN_ID}" \
        "${SHUFFLE_SEED}" \
        >> "${CAMPAIGN_DIR}/shuffle-seeds.csv"

    echo
    echo "============================================================"
    echo " Run ${RUN_ID}/${RUNS}"
    echo " Shuffle seed: ${SHUFFLE_SEED}"
    echo " Output      : ${OUT}"
    echo "============================================================"

    CXX=g++ ./benchmark/scripts/run-benchmark.sh \
        --cpu 17 \
        --repetitions 200 \
        --warmup 50 \
        --batch 1 \
        --pmu-profile none \
        --degrees 2048,4096,8192,16384,32768 \
        --bits 20,25,30,40,50,52,54,56,58,60,61,62 \
        --kernels all \
        --axhel-cpu native \
        --shuffle-scenarios \
        --shuffle-seed "${SHUFFLE_SEED}" \
        --output "${OUT}" \
        2>&1 | tee "${CAMPAIGN_DIR}/run-${RUN_ID}.log"

    sync
    sleep 2
done

echo
echo "============================================================"
echo " Campaign completed"
echo " Output: ${CAMPAIGN_DIR}"
echo "============================================================"
