// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include "axhel/util/compiler.hpp"

#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>

namespace unipi {
namespace axhel {

    void EltwiseAddModSVE(uint64_t *result, const uint64_t *operand1, const uint64_t *operand2, uint64_t n, uint64_t modulus) {
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

            const svuint64_t sum0 = svadd_u64_x(pg_all, a0, b0);
            const svuint64_t sum1 = svadd_u64_x(pg_all, a1, b1);
            const svuint64_t sum2 = svadd_u64_x(pg_all, a2, b2);
            const svuint64_t sum3 = svadd_u64_x(pg_all, a3, b3);

            const svuint64_t reduced0 = svsub_u64_x(pg_all, sum0, modulus_vec);
            const svuint64_t reduced1 = svsub_u64_x(pg_all, sum1, modulus_vec);
            const svuint64_t reduced2 = svsub_u64_x(pg_all, sum2, modulus_vec);
            const svuint64_t reduced3 = svsub_u64_x(pg_all, sum3, modulus_vec);

            // Operands are assumed to be reduced modulo q, so sum is in [0, 2q).
            // Branchless conditional subtraction: min(sum, sum - modulus).
            const svuint64_t r0 = svmin_u64_x(pg_all, sum0, reduced0);
            const svuint64_t r1 = svmin_u64_x(pg_all, sum1, reduced1);
            const svuint64_t r2 = svmin_u64_x(pg_all, sum2, reduced2);
            const svuint64_t r3 = svmin_u64_x(pg_all, sum3, reduced3);

            svst1_u64(pg_all, result + i, r0);
            svst1_u64(pg_all, result + i + lanes, r1);
            svst1_u64(pg_all, result + i + 2 * lanes, r2);
            svst1_u64(pg_all, result + i + 3 * lanes, r3);
        }

        // Handle the remaining elements with predication.
        while (i < n) {
            const svbool_t pg = svwhilelt_b64(i, n);

            const svuint64_t a = svld1_u64(pg, operand1 + i);
            const svuint64_t b = svld1_u64(pg, operand2 + i);
            
            const svuint64_t sum = svadd_u64_x(pg, a, b);
            
            const svuint64_t reduced = svsub_u64_x(pg, sum, modulus_vec);
            
            const svuint64_t r = svmin_u64_x(pg, sum, reduced);

            svst1_u64(pg, result + i, r);

            i += lanes;
        }
    }


    void EltwiseAddModSVE(uint64_t *result, const uint64_t *operand1, uint64_t operand2, uint64_t n, uint64_t modulus) {
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

            const svuint64_t sum0 = svadd_u64_x(pg_all, a0, operand2_vec);
            const svuint64_t sum1 = svadd_u64_x(pg_all, a1, operand2_vec);
            const svuint64_t sum2 = svadd_u64_x(pg_all, a2, operand2_vec);
            const svuint64_t sum3 = svadd_u64_x(pg_all, a3, operand2_vec);

            const svuint64_t reduced0 = svsub_u64_x(pg_all, sum0, modulus_vec);
            const svuint64_t reduced1 = svsub_u64_x(pg_all, sum1, modulus_vec);
            const svuint64_t reduced2 = svsub_u64_x(pg_all, sum2, modulus_vec);
            const svuint64_t reduced3 = svsub_u64_x(pg_all, sum3, modulus_vec);

            // Operands are assumed to be reduced modulo q, so sum is in [0, 2q).
            // Branchless conditional subtraction: min(sum, sum - modulus).
            const svuint64_t r0 = svmin_u64_x(pg_all, sum0, reduced0);
            const svuint64_t r1 = svmin_u64_x(pg_all, sum1, reduced1);
            const svuint64_t r2 = svmin_u64_x(pg_all, sum2, reduced2);
            const svuint64_t r3 = svmin_u64_x(pg_all, sum3, reduced3);

            svst1_u64(pg_all, result + i, r0);
            svst1_u64(pg_all, result + i + lanes, r1);
            svst1_u64(pg_all, result + i + 2 * lanes, r2);
            svst1_u64(pg_all, result + i + 3 * lanes, r3);
        }

        // Handle the remaining elements with predication.
        while (i < n) {
            const svbool_t pg = svwhilelt_b64(i, n);

            const svuint64_t a = svld1_u64(pg, operand1 + i);

            const svuint64_t sum = svadd_u64_x(pg, a, operand2_vec);

            const svuint64_t reduced = svsub_u64_x(pg, sum, modulus_vec);

            const svuint64_t r = svmin_u64_x(pg, sum, reduced);

            svst1_u64(pg, result + i, r);

            i += lanes;
        }
    }


}
}

#endif
