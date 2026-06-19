// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#if defined(AXHEL_HAS_SVE) || defined(AXHEL_HAS_SVE2)

#include <arm_sve.h>
#include <stdint.h>

namespace unipi {
namespace axhel {

    //@brief Compile-time 128-bit right shift producing the lower 64-bit result
    template <uint64_t Shift>
    inline svuint64_t ShiftRight128LowPart(svbool_t pg, svuint64_t hi, svuint64_t lo) {
    if constexpr (Shift == 0) {
            return lo;
        } else {
            svuint64_t lo_part = svlsr_n_u64_x(pg, lo, Shift);
            svuint64_t hi_part = svlsl_n_u64_x(pg, hi, 64 - Shift);
            return svorr_u64_x(pg, lo_part, hi_part);
        }
    }

    // @brief SVE-only high 64 bits of a 64x64 unsigned multiplication.
    inline svuint64_t MulHighU64SVE(svbool_t pg, svuint64_t x, svuint64_t y) {
        const svuint64_t mask32 = svdup_n_u64(0xFFFFFFFFULL);

        svuint64_t x_lo = svand_u64_x(pg, x, mask32);
        svuint64_t x_hi = svlsr_n_u64_x(pg, x, 32);

        svuint64_t y_lo = svand_u64_x(pg, y, mask32);
        svuint64_t y_hi = svlsr_n_u64_x(pg, y, 32);

        svuint64_t p0 = svmul_u64_x(pg, x_lo, y_lo);
        svuint64_t p1 = svmul_u64_x(pg, x_lo, y_hi);
        svuint64_t p2 = svmul_u64_x(pg, x_hi, y_lo);
        svuint64_t p3 = svmul_u64_x(pg, x_hi, y_hi);

        // Carry from the low 64-bit half.
        svuint64_t middle = svadd_u64_x(pg, svlsr_n_u64_x(pg, p0, 32), svand_u64_x(pg, p1, mask32));
        middle = svadd_u64_x(pg, middle, svand_u64_x(pg, p2, mask32));
        svuint64_t carry = svlsr_n_u64_x(pg, middle, 32);

        // High 64 bits of x*y.
        svuint64_t hi =svadd_u64_x(pg, p3, svlsr_n_u64_x(pg, p1, 32));
        hi = svadd_u64_x(pg, hi, svlsr_n_u64_x(pg, p2, 32));
        hi =  svadd_u64_x(pg, hi, carry);

        return hi;
    }

    // @brief SVE-only 64x64 -> 128 unsigned multiplication.
    inline void MulU64ToU128SVE(svbool_t pg, svuint64_t x, svuint64_t y, svuint64_t* hi, svuint64_t* lo) {
        *lo = svmul_u64_x(pg, x, y);
        *hi = MulHighU64SVE(pg, x, y);
    }


}
}
#endif