#!/usr/bin/env bash

# Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

RUNS="${1:-30}"
BASE_SHUFFLE_SEED="${2:-2026082900}"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT_DIR}"

TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
CAMPAIGN_DIR="benchmark/results/${TIMESTAMP}-all-kernels-definitive"

ELEMENTWISE_KERNELS="add_vv,add_vs,sub_vv,sub_vs,mul_f1,mul_f2,mul_f4,fma_f1,fma_f2,fma_f4,fma_f8,reduce_2to1,reduce_4to2,reduce_4to1,reduce_barrett"
NTT_KERNELS="ntt_forward,ntt_forward_lazy,ntt_inverse,ntt_inverse_lazy"

mkdir -p "${CAMPAIGN_DIR}"

echo "AXHEL all-kernels randomized optimization campaign"
echo "Runs              : ${RUNS}"
echo "Base shuffle seed : ${BASE_SHUFFLE_SEED}"
echo "CPU               : 17"
echo "Repetitions       : 200"
echo "Warmup            : 50"
echo "Element-wise batch: 16"
echo "NTT batch         : 1"
echo "PMU profile       : none"
echo "Degrees           : 2048,4096,8192,16384,32768,65536"
echo "qi bits           : 20,25,30,40,50,52,54,56,58,60,61,62"
echo "AXHEL CPU target   : native"
echo "Scenario order    : randomized"
echo "Output            : ${CAMPAIGN_DIR}"
echo "Git commit        : $(git rev-parse HEAD)"
echo

# Snapshot exact source state used for the campaign.
git rev-parse HEAD > "${CAMPAIGN_DIR}/git-commit.txt"
git status --short > "${CAMPAIGN_DIR}/git-status.txt"
git diff > "${CAMPAIGN_DIR}/git-diff.patch"

cat > "${CAMPAIGN_DIR}/campaign-config.txt" <<EOF
runs=${RUNS}
base_shuffle_seed=${BASE_SHUFFLE_SEED}
cpu=17
repetitions=200
warmup=50
elementwise_batch=16
ntt_batch=1
pmu_profile=none
degrees=2048,4096,8192,16384,32768,65536
bits=20,25,30,40,50,52,54,56,58,60,61,62
elementwise_kernels=${ELEMENTWISE_KERNELS}
ntt_kernels=${NTT_KERNELS}
axhel_cpu=native
scenario_order=randomized
raw_samples=enabled
EOF

printf "run,shuffle_seed\n" > "${CAMPAIGN_DIR}/shuffle-seeds.csv"

for RUN in $(seq 1 "${RUNS}"); do
    RUN_ID="$(printf "%02d" "${RUN}")"
    SHUFFLE_SEED=$((BASE_SHUFFLE_SEED + RUN))
    RUN_DIR="${CAMPAIGN_DIR}/run-${RUN_ID}"

    mkdir -p "${RUN_DIR}"

    printf "%s,%s\n" \
        "${RUN_ID}" \
        "${SHUFFLE_SEED}" \
        >> "${CAMPAIGN_DIR}/shuffle-seeds.csv"

    echo
    echo "============================================================"
    echo " Run ${RUN_ID}/${RUNS}"
    echo " Shuffle seed: ${SHUFFLE_SEED}"
    echo "============================================================"

    echo
    echo "---- Element-wise / reductions : batch=16 ----"
    CXX=g++ ./benchmark/scripts/run-benchmark.sh \
        --cpu 17 \
        --repetitions 200 \
        --warmup 50 \
        --batch 16 \
        --pmu-profile none \
        --degrees 2048,4096,8192,16384,32768,65536 \
        --bits 20,25,30,40,50,52,54,56,58,60,61,62 \
        --kernels "${ELEMENTWISE_KERNELS}" \
        --axhel-cpu native \
        --shuffle-scenarios \
        --shuffle-seed "${SHUFFLE_SEED}" \
        --output "${RUN_DIR}/elementwise" \
        2>&1 | tee "${CAMPAIGN_DIR}/run-${RUN_ID}-elementwise.log"

    sync

    echo
    echo "---- NTT / INTT : batch=1 ----"
    CXX=g++ ./benchmark/scripts/run-benchmark.sh \
        --cpu 17 \
        --repetitions 200 \
        --warmup 50 \
        --batch 1 \
        --pmu-profile none \
        --degrees 2048,4096,8192,16384,32768,65536 \
        --bits 20,25,30,40,50,52,54,56,58,60,61,62 \
        --kernels "${NTT_KERNELS}" \
        --axhel-cpu native \
        --shuffle-scenarios \
        --shuffle-seed "${SHUFFLE_SEED}" \
        --output "${RUN_DIR}/ntt" \
        2>&1 | tee "${CAMPAIGN_DIR}/run-${RUN_ID}-ntt.log"

    sync
    sleep 2
done

echo
echo "============================================================"
echo " Campaign completed"
echo " Output: ${CAMPAIGN_DIR}"
echo "============================================================"
