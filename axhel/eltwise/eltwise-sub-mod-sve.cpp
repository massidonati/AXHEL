// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include "axhel/util/compiler.hpp"

#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>

namespace unipi {
namespace axhel {

    void EltwiseSubModSVE(uint64_t *result, const uint64_t *operand1, const uint64_t *operand2, uint64_t n, uint64_t modulus) {
        const uint64_t lanes = static_cast<uint64_t>(svcntd());
        const uint64_t block = 4 * lanes;

        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);

        uint64_t i = 0;

        // Main loop: process four independent full-width SVE vectors.
        for (; i + block <= n; i += block) {
            const svuint64_t a0 = svld1_u64(pg_all, operand1 + i);
            const svuint64_t a1 = svld1_u64(pg_all, operand1 + i + lanes);
            const svuint64_t a2 = svld1_u64(pg_all, operand1 + i + 2 * lanes);
            const svuint64_t a3 = svld1_u64(pg_all, operand1 + i + 3 * lanes);

            const svuint64_t b0 = svld1_u64(pg_all, operand2 + i);
            const svuint64_t b1 = svld1_u64(pg_all, operand2 + i + lanes);
            const svuint64_t b2 = svld1_u64(pg_all, operand2 + i + 2 * lanes);
            const svuint64_t b3 = svld1_u64(pg_all, operand2 + i + 3 * lanes);

            const svuint64_t d0 = svsub_u64_x(pg_all, a0, b0);
            const svuint64_t d1 = svsub_u64_x(pg_all, a1, b1);
            const svuint64_t d2 = svsub_u64_x(pg_all, a2, b2);
            const svuint64_t d3 = svsub_u64_x(pg_all, a3, b3);

            const svuint64_t c0 = svadd_u64_x(pg_all, d0, modulus_vec);
            const svuint64_t c1 = svadd_u64_x(pg_all, d1, modulus_vec);
            const svuint64_t c2 = svadd_u64_x(pg_all, d2, modulus_vec);
            const svuint64_t c3 = svadd_u64_x(pg_all, d3, modulus_vec);

            // Operands are assumed to be reduced modulo q, so the mathematical difference is in (-q, q).
            // Branchless modular correction: min(diff, diff + modulus).
            const svuint64_t r0 = svmin_u64_x(pg_all, d0, c0);
            const svuint64_t r1 = svmin_u64_x(pg_all, d1, c1);
            const svuint64_t r2 = svmin_u64_x(pg_all, d2, c2);
            const svuint64_t r3 = svmin_u64_x(pg_all, d3, c3);

            svst1_u64(pg_all, result + i,             r0);
            svst1_u64(pg_all, result + i + lanes,     r1);
            svst1_u64(pg_all, result + i + 2 * lanes, r2);
            svst1_u64(pg_all, result + i + 3 * lanes, r3);
        }

        // Handle the remaining elements with predication.
        for (; i < n; i += lanes) {
            const svbool_t pg = svwhilelt_b64(i, n);

            const svuint64_t a = svld1_u64(pg, operand1 + i);
            const svuint64_t b = svld1_u64(pg, operand2 + i);

            const svuint64_t d = svsub_u64_x(pg, a, b);

            const svuint64_t corrected = svadd_u64_x(pg, d, modulus_vec);

            const svuint64_t r = svmin_u64_x(pg, d, corrected);

            svst1_u64(pg, result + i, r);
        }
    }


    void EltwiseSubModSVE(uint64_t *result, const uint64_t *operand1, uint64_t operand2, uint64_t n, uint64_t modulus) {
        const uint64_t lanes = static_cast<uint64_t>(svcntd());
        const uint64_t block = 4 * lanes;

        const svbool_t pg_all = svptrue_b64();
        const svuint64_t operand2_vec = svdup_n_u64(operand2);
        const svuint64_t modulus_vec = svdup_n_u64(modulus);

        uint64_t i = 0;

        // Main loop: process four independent full-width SVE vectors.
        for (; i + block <= n; i += block) {
            const svuint64_t a0 = svld1_u64(pg_all, operand1 + i);
            const svuint64_t a1 = svld1_u64(pg_all, operand1 + i + lanes);
            const svuint64_t a2 = svld1_u64(pg_all, operand1 + i + 2 * lanes);
            const svuint64_t a3 = svld1_u64(pg_all, operand1 + i + 3 * lanes);

            const svuint64_t d0 = svsub_u64_x(pg_all, a0, operand2_vec);
            const svuint64_t d1 = svsub_u64_x(pg_all, a1, operand2_vec);
            const svuint64_t d2 = svsub_u64_x(pg_all, a2, operand2_vec);
            const svuint64_t d3 = svsub_u64_x(pg_all, a3, operand2_vec);

            const svuint64_t c0 = svadd_u64_x(pg_all, d0, modulus_vec);
            const svuint64_t c1 = svadd_u64_x(pg_all, d1, modulus_vec);
            const svuint64_t c2 = svadd_u64_x(pg_all, d2, modulus_vec);
            const svuint64_t c3 = svadd_u64_x(pg_all, d3, modulus_vec);

            // Operands are assumed to be reduced modulo q, so the mathematical difference is in (-q, q).
            // Branchless modular correction: min(diff, diff + modulus).
            const svuint64_t r0 = svmin_u64_x(pg_all, d0, c0);
            const svuint64_t r1 = svmin_u64_x(pg_all, d1, c1);
            const svuint64_t r2 = svmin_u64_x(pg_all, d2, c2);
            const svuint64_t r3 = svmin_u64_x(pg_all, d3, c3);

            svst1_u64(pg_all, result + i,             r0);
            svst1_u64(pg_all, result + i + lanes,     r1);
            svst1_u64(pg_all, result + i + 2 * lanes, r2);
            svst1_u64(pg_all, result + i + 3 * lanes, r3);
        }

        // Handle the remaining elements with predication.
        for (; i < n; i += lanes) {
            const svbool_t pg = svwhilelt_b64(i, n);

            const svuint64_t a = svld1_u64(pg, operand1 + i);

            const svuint64_t d = svsub_u64_x(pg, a, operand2_vec);

            const svuint64_t corrected = svadd_u64_x(pg, d, modulus_vec);

            const svuint64_t r = svmin_u64_x(pg, d, corrected);

            svst1_u64(pg, result + i, r);
        }
    }

    
}
}

#endif
