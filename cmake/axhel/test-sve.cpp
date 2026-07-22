// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <arm_sve.h>
#include <cstdint>

#ifndef __ARM_FEATURE_SVE
#error "SVE not enabled for this target"
#endif

int main(void) {
    const svbool_t pg = svptrue_b64();
    const svuint64_t one = svdup_u64(1);
    const svuint64_t two = svdup_u64(2);
    const svuint64_t sum = svadd_u64_x(pg, one, two);
    const svuint64_t high = svmulh_u64_x(pg, one, two);
    const std::uint64_t sum_result = svlastb_u64(pg, sum);
    const std::uint64_t high_result = svlastb_u64(pg, high);

    return (sum_result == 3 && high_result == 0) ? 0 : 1;
}