// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include "eltwise/eltwise-reduce-mod-sve.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/sve-arith.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "axhel/util/compiler.hpp"


#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>

namespace unipi {
namespace axhel {


    // SVE Barrett reduction using a compile-time shift.
    template <uint64_t Shift>
    inline void BarrettReduceUInt64SVEKernel(uint64_t *res, const uint64_t *op, uint64_t n, uint64_t mod, uint64_t barr_factor) {
        const uint64_t lanes = svcntd();
        const uint64_t block = 4 * lanes;

        const svbool_t pg_all = svptrue_b64();
        const svuint64_t vmod = svdup_n_u64(mod);
        const svuint64_t v2mod = svdup_n_u64(2 * mod);

        // Required by the ReduceInputSVE<4>() interface.
        const svuint64_t v4mod = svdup_n_u64(4 * mod);
        const svuint64_t vbarr = svdup_n_u64(barr_factor);

        uint64_t i = 0;

        // Main loop: process four independent full-width SVE vectors.
        for (; n - i >= block; i += block) {
            const svuint64_t op0 = svld1_u64(pg_all, op + i);
            const svuint64_t op1 = svld1_u64(pg_all, op + i + lanes);
            const svuint64_t op2 = svld1_u64(pg_all, op + i + 2 * lanes);
            const svuint64_t op3 = svld1_u64(pg_all, op + i + 3 * lanes);

            // Barrett quotient approximation: c = op >> Shift.
            const svuint64_t c0 = svlsr_n_u64_x(pg_all, op0, Shift);
            const svuint64_t c1 = svlsr_n_u64_x(pg_all, op1, Shift);
            const svuint64_t c2 = svlsr_n_u64_x(pg_all, op2, Shift);
            const svuint64_t c3 = svlsr_n_u64_x(pg_all, op3, Shift);

            // Approximate quotient: q_hat = high64(c * barr_factor).
            const svuint64_t q0 = svmulh_u64_x(pg_all, c0, vbarr);
            const svuint64_t q1 = svmulh_u64_x(pg_all, c1, vbarr);
            const svuint64_t q2 = svmulh_u64_x(pg_all, c2, vbarr);
            const svuint64_t q3 = svmulh_u64_x(pg_all, c3, vbarr);

            const svuint64_t qm0 = svmul_u64_x(pg_all, q0, vmod);
            const svuint64_t qm1 = svmul_u64_x(pg_all, q1, vmod);
            const svuint64_t qm2 = svmul_u64_x(pg_all, q2, vmod);
            const svuint64_t qm3 = svmul_u64_x(pg_all, q3, vmod);

            // Barrett residual.
            svuint64_t z0 = svsub_u64_x(pg_all, op0, qm0);
            svuint64_t z1 = svsub_u64_x(pg_all, op1, qm1);
            svuint64_t z2 = svsub_u64_x(pg_all, op2, qm2);
            svuint64_t z3 = svsub_u64_x(pg_all, op3, qm3);

            // Final correction: [0, 4q) -> [0, q).
            z0 = ReduceInputSVE<4>(pg_all, z0, vmod, v2mod, v4mod);
            z1 = ReduceInputSVE<4>(pg_all, z1, vmod, v2mod, v4mod);
            z2 = ReduceInputSVE<4>(pg_all, z2, vmod, v2mod, v4mod);
            z3 = ReduceInputSVE<4>(pg_all, z3, vmod, v2mod, v4mod);

            svst1_u64(pg_all, res + i, z0);
            svst1_u64(pg_all, res + i + lanes, z1);
            svst1_u64(pg_all, res + i + 2 * lanes, z2);
            svst1_u64(pg_all, res + i + 3 * lanes, z3);    
        }

        // Handle the remaining elements with predication.
        while (i < n) {
            const svbool_t pg =svwhilelt_b64(i, n);

            const svuint64_t op0 = svld1_u64(pg, op + i);

            const svuint64_t c0 = svlsr_n_u64_x( pg, op0, Shift);

            const svuint64_t q0 = svmulh_u64_x(pg, c0, vbarr);

            const svuint64_t qm0 = svmul_u64_x(pg, q0, vmod);

            svuint64_t z0 = svsub_u64_x(pg, op0, qm0);

            z0 = ReduceInputSVE<4>(pg, z0, vmod, v2mod, v4mod);

            svst1_u64(pg, res + i, z0);

            i += lanes;
        }
    }


    // Dispatches SVE Barrett reduction to the kernel matching the runtime shift value.
    void DispatchBarrettReduceUInt64SVEKernel(uint64_t *res, const uint64_t *op, uint64_t n, uint64_t mod, uint64_t barr_factor, uint64_t prod_right_shift) {
    #define AXHEL_DISPATCH_SHIFT(SHIFT_VALUE)              \
        case SHIFT_VALUE:                                  \
            BarrettReduceUInt64SVEKernel<SHIFT_VALUE>(     \
                res,                                       \
                op,                                        \
                n,                                         \
                mod,                                       \
                barr_factor);                              \
            break

        switch (prod_right_shift) {
            AXHEL_DISPATCH_SHIFT(1);
            AXHEL_DISPATCH_SHIFT(2);
            AXHEL_DISPATCH_SHIFT(3);
            AXHEL_DISPATCH_SHIFT(4);
            AXHEL_DISPATCH_SHIFT(5);
            AXHEL_DISPATCH_SHIFT(6);
            AXHEL_DISPATCH_SHIFT(7);
            AXHEL_DISPATCH_SHIFT(8);
            AXHEL_DISPATCH_SHIFT(9);
            AXHEL_DISPATCH_SHIFT(10);
            AXHEL_DISPATCH_SHIFT(11);
            AXHEL_DISPATCH_SHIFT(12);
            AXHEL_DISPATCH_SHIFT(13);
            AXHEL_DISPATCH_SHIFT(14);
            AXHEL_DISPATCH_SHIFT(15);
            AXHEL_DISPATCH_SHIFT(16);
            AXHEL_DISPATCH_SHIFT(17);
            AXHEL_DISPATCH_SHIFT(18);
            AXHEL_DISPATCH_SHIFT(19);
            AXHEL_DISPATCH_SHIFT(20);
            AXHEL_DISPATCH_SHIFT(21);
            AXHEL_DISPATCH_SHIFT(22);
            AXHEL_DISPATCH_SHIFT(23);
            AXHEL_DISPATCH_SHIFT(24);
            AXHEL_DISPATCH_SHIFT(25);
            AXHEL_DISPATCH_SHIFT(26);
            AXHEL_DISPATCH_SHIFT(27);
            AXHEL_DISPATCH_SHIFT(28);
            AXHEL_DISPATCH_SHIFT(29);
            AXHEL_DISPATCH_SHIFT(30);
            AXHEL_DISPATCH_SHIFT(31);
            AXHEL_DISPATCH_SHIFT(32);
            AXHEL_DISPATCH_SHIFT(33);
            AXHEL_DISPATCH_SHIFT(34);
            AXHEL_DISPATCH_SHIFT(35);
            AXHEL_DISPATCH_SHIFT(36);
            AXHEL_DISPATCH_SHIFT(37);
            AXHEL_DISPATCH_SHIFT(38);
            AXHEL_DISPATCH_SHIFT(39);
            AXHEL_DISPATCH_SHIFT(40);
            AXHEL_DISPATCH_SHIFT(41);
            AXHEL_DISPATCH_SHIFT(42);
            AXHEL_DISPATCH_SHIFT(43);
            AXHEL_DISPATCH_SHIFT(44);
            AXHEL_DISPATCH_SHIFT(45);
            AXHEL_DISPATCH_SHIFT(46);
            AXHEL_DISPATCH_SHIFT(47);
            AXHEL_DISPATCH_SHIFT(48);
            AXHEL_DISPATCH_SHIFT(49);
            AXHEL_DISPATCH_SHIFT(50);
            AXHEL_DISPATCH_SHIFT(51);
            AXHEL_DISPATCH_SHIFT(52);
            AXHEL_DISPATCH_SHIFT(53);
            AXHEL_DISPATCH_SHIFT(54);
            AXHEL_DISPATCH_SHIFT(55);
            AXHEL_DISPATCH_SHIFT(56);
            AXHEL_DISPATCH_SHIFT(57);
            AXHEL_DISPATCH_SHIFT(58);
            AXHEL_DISPATCH_SHIFT(59);
            AXHEL_DISPATCH_SHIFT(60);

            default:
                ; // TODO error
        }

    #undef AXHEL_DISPATCH_SHIFT
    }



    namespace {

        /*
        * Bounded SVE reductions:
        *
        *     <2, 1>: [0, 2q) -> [0, q)
        *     <4, 2>: [0, 4q) -> [0, 2q)
        *     <4, 1>: [0, 4q) -> [0, q)
        *
        * In and Out are compile-time bounding factors.
        */
        template <int In, int Out>
        void EltwiseReduceModBoundedSVE(uint64_t *res, const uint64_t *op, uint64_t n, uint64_t mod) {

            const uint64_t lanes = svcntd();
            const uint64_t block = 4 * lanes;

            const svbool_t pg_all = svptrue_b64();
            const svuint64_t vmod = svdup_n_u64(mod);
            const svuint64_t v2mod = svdup_n_u64(2 * mod);

            uint64_t i = 0;

            // Main loop: process four independent full-width SVE vectors.
            for (; n - i >= block; i += block) {
                svuint64_t x0 = svld1_u64(pg_all, op + i);
                svuint64_t x1 = svld1_u64(pg_all, op + i + lanes);
                svuint64_t x2 = svld1_u64(pg_all, op + i + 2 * lanes);
                svuint64_t x3 = svld1_u64(pg_all, op + i + 3 * lanes);

                // [0, 4q) -> [0, 2q).
                if constexpr (In == 4) {
                    x0 = ReduceModFactor4To2SVE(pg_all, x0, v2mod);
                    x1 = ReduceModFactor4To2SVE(pg_all, x1, v2mod);
                    x2 = ReduceModFactor4To2SVE(pg_all, x2, v2mod);
                    x3 = ReduceModFactor4To2SVE(pg_all, x3, v2mod);
                }

                // [0, 2q) -> [0, q).
                if constexpr (Out == 1) {
                    x0 = ReduceModFactor2To1SVE(pg_all, x0, vmod);
                    x1 = ReduceModFactor2To1SVE(pg_all, x1, vmod);
                    x2 = ReduceModFactor2To1SVE(pg_all, x2, vmod);
                    x3 = ReduceModFactor2To1SVE(pg_all, x3, vmod);
                }

                svst1_u64(pg_all, res + i, x0);
                svst1_u64(pg_all, res + i + lanes, x1);
                svst1_u64(pg_all, res + i + 2 * lanes, x2);
                svst1_u64(pg_all, res + i + 3 * lanes, x3);
            }

            /// Handle the remaining elements with predication.
            while (i < n) {
                const svbool_t pg = svwhilelt_b64(i, n);

                svuint64_t x = svld1_u64(pg, op + i);

                if constexpr (In == 4) {
                    x = ReduceModFactor4To2SVE(pg, x, v2mod);
                }

                if constexpr (Out == 1) {
                    x = ReduceModFactor2To1SVE(pg, x, vmod);
                }

                svst1_u64(pg, res + i, x);

                i += lanes;
            }
        }

    } // namespace



    void EltwiseReduceModSVE(uint64_t *res, const uint64_t *op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor) {
        const uint64_t lanes = svcntd();

        // Identity reductions.
        if (in_mod_factor == 1 || (in_mod_factor == 2 && out_mod_factor == 2)) {
            // No copy is required for in-place operation.
            if (res == op) {
                return;
            }

            // Predicated copy.
            for (uint64_t i = 0; i < n; i += lanes) {
                const svbool_t pg = svwhilelt_b64(i, n);

                const svuint64_t x = svld1_u64(pg, op + i);

                svst1_u64(pg, res + i, x);
            }

            return;
        }

        // [0, 2q) -> [0, q).
        if (in_mod_factor == 2 && out_mod_factor == 1) {
            EltwiseReduceModBoundedSVE<2, 1>(res, op, n, mod);
            return;
        }

        // [0, 4q) -> [0, q).
        if (in_mod_factor == 4 && out_mod_factor == 1) {
            EltwiseReduceModBoundedSVE<4, 1>(res, op, n, mod);
            return;
        }

        // [0, 4q) -> [0, 2q)
        if (in_mod_factor == 4 && out_mod_factor == 2) {
            EltwiseReduceModBoundedSVE<4, 2>(res, op, n, mod);
            return;
        }

        // Barrett reduction for the general input range [0, mod^2).
        if (in_mod_factor == mod) {
            constexpr int64_t beta = -2;
            constexpr int64_t alpha = 62;

            const uint64_t ceil_log_mod = BitWidth(mod);

            const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

            const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

            DispatchBarrettReduceUInt64SVEKernel(res, op, n, mod, barr_factor, prod_right_shift);

            return;
        }
    }


}
}

#endif