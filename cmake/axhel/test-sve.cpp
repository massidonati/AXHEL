// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <arm_sve.h>
#include <cstdint>

#ifndef __ARM_FEATURE_SVE
#error "SVE not enabled for this target"
#endif

int main(void) {
    const svbool_t pg = svptrue_b64();
    const svint64_t one = svdup_s64(1);
    const svint64_t two = svdup_s64(2);
    const svint64_t sum = svadd_s64_x(pg, one, two);
    const std::int64_t result = svlastb_s64(pg, sum);
    constexpr std::int64_t expected = 3;
    return (result == expected) ? 0 : 1;
}