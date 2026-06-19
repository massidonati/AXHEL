// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <arm_sve.h>

#ifndef __ARM_FEATURE_SVE2
#error "SVE2 not enabled for this target"
#endif

int main(void) {
    svint32_t a = svdup_s32(1000);
    svint32_t b = svdup_s32(2000);
    // SVE2 intrinsic
    svint32_t r = svqdmulh_s32(a, b);
    return svlastb_s32(svptrue_b32(), r) >= 0 ? 0 : 1;
}