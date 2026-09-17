// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include "eltwise/eltwise-fma-mod-sve.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/sve-arith.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "axhel/util/compiler.hpp"

#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>

namespace unipi {
namespace axhel {


    // Performs SVE element-wise modular multiply-add using a compile-time shift.
    template <int ModFactor, uint64_t Shift>
    void EltwiseFMAModSVEKernel(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t barr_factor) {
        const uint64_t lanes = svcntd();
        const uint64_t block = 4 * lanes;

        const svbool_t pg_all = svptrue_b64();
        const svuint64_t vmod = svdup_n_u64(mod);
        const svuint64_t v2mod = svdup_n_u64(2 * mod);
        const svuint64_t v4mod = svdup_n_u64(4 * mod);
        const svuint64_t vbarr = svdup_n_u64(barr_factor);

        // Reduce scalar operand to [0, q).
        const uint64_t op2_red = ReduceInputNative<ModFactor>(op2, mod);
        const svuint64_t yop = svdup_n_u64(op2_red);

        uint64_t i = 0;

        // Main loop: process four independent full-width SVE vectors.
        for (; n - i >= block; i += block) {

            svuint64_t x0 = svld1_u64(pg_all, op1 + i);
            svuint64_t x1 = svld1_u64(pg_all, op1 + i + lanes);
            svuint64_t x2 = svld1_u64(pg_all, op1 + i + 2 * lanes);
            svuint64_t x3 = svld1_u64(pg_all, op1 + i + 3 * lanes);

            // Reduce inputs to [0, q).
            if constexpr (ModFactor != 1) {
                x0 = ReduceInputSVE<ModFactor>(pg_all, x0, vmod, v2mod, v4mod);
                x1 = ReduceInputSVE<ModFactor>(pg_all, x1, vmod, v2mod, v4mod);
                x2 = ReduceInputSVE<ModFactor>(pg_all, x2, vmod, v2mod, v4mod);
                x3 = ReduceInputSVE<ModFactor>(pg_all, x3, vmod, v2mod, v4mod);
            }

            // Full 64x64 -> 128-bit products.
            svuint64_t prod_hi0;
            svuint64_t prod_lo0;

            svuint64_t prod_hi1;
            svuint64_t prod_lo1;

            svuint64_t prod_hi2;
            svuint64_t prod_lo2;

            svuint64_t prod_hi3;
            svuint64_t prod_lo3;

            MulU64ToU128SVE(pg_all, x0, yop, &prod_hi0, &prod_lo0);
            MulU64ToU128SVE(pg_all, x1, yop, &prod_hi1, &prod_lo1);
            MulU64ToU128SVE(pg_all, x2, yop, &prod_hi2, &prod_lo2);
            MulU64ToU128SVE(pg_all, x3, yop, &prod_hi3, &prod_lo3);

            // Extract the low 64 bits of the shifted 128-bit products.
            const svuint64_t c0 = ShiftRight128LowPart<Shift>(pg_all, prod_hi0, prod_lo0);
            const svuint64_t c1 = ShiftRight128LowPart<Shift>(pg_all, prod_hi1, prod_lo1);
            const svuint64_t c2 = ShiftRight128LowPart<Shift>(pg_all, prod_hi2, prod_lo2);
            const svuint64_t c3 = ShiftRight128LowPart<Shift>(pg_all, prod_hi3, prod_lo3);

            // Barrett quotient approximations.
            const svuint64_t q_hat0 = svmulh_u64_x(pg_all, c0, vbarr);
            const svuint64_t q_hat1 = svmulh_u64_x(pg_all, c1, vbarr);
            const svuint64_t q_hat2 = svmulh_u64_x(pg_all, c2, vbarr);
            const svuint64_t q_hat3 = svmulh_u64_x(pg_all, c3, vbarr);

            // Barrett residuals: z = prod_lo - q_hat * mod.
            const svuint64_t q_mul0 = svmul_u64_x(pg_all, q_hat0, vmod);
            const svuint64_t q_mul1 = svmul_u64_x(pg_all, q_hat1, vmod);
            const svuint64_t q_mul2 = svmul_u64_x(pg_all, q_hat2, vmod);
            const svuint64_t q_mul3 = svmul_u64_x(pg_all, q_hat3, vmod);

            svuint64_t z0 = svsub_u64_x(pg_all, prod_lo0, q_mul0);
            svuint64_t z1 = svsub_u64_x(pg_all, prod_lo1, q_mul1);
            svuint64_t z2 = svsub_u64_x(pg_all, prod_lo2, q_mul2);
            svuint64_t z3 = svsub_u64_x(pg_all, prod_lo3, q_mul3);

            // Final correction: [0, 4q) -> [0, q).
            z0 = ReduceInputSVE<4>(pg_all, z0, vmod, v2mod, v4mod);
            z1 = ReduceInputSVE<4>(pg_all, z1, vmod, v2mod, v4mod);
            z2 = ReduceInputSVE<4>(pg_all, z2, vmod, v2mod, v4mod);
            z3 = ReduceInputSVE<4>(pg_all, z3, vmod, v2mod, v4mod);

            // Add op3 when provided.
            if (op3 != nullptr) {

                svuint64_t add0 = svld1_u64(pg_all, op3 + i);
                svuint64_t add1 = svld1_u64(pg_all, op3 + i + lanes);
                svuint64_t add2 = svld1_u64(pg_all, op3 + i + 2 * lanes);
                svuint64_t add3 = svld1_u64(pg_all, op3 + i + 3 * lanes);

                // Reduce addends to [0, q).
                if constexpr (ModFactor != 1) {
                    add0 = ReduceInputSVE<ModFactor>(pg_all, add0, vmod, v2mod, v4mod);
                    add1 = ReduceInputSVE<ModFactor>(pg_all, add1, vmod, v2mod, v4mod);
                    add2 = ReduceInputSVE<ModFactor>(pg_all, add2, vmod, v2mod, v4mod);
                    add3 = ReduceInputSVE<ModFactor>(pg_all, add3, vmod, v2mod, v4mod);
                }

                // Modular addition: [0, q) + [0, q) -> [0, 2q) -> [0, q).
                z0 = svadd_u64_x(pg_all, z0, add0);
                z1 = svadd_u64_x(pg_all, z1, add1);
                z2 = svadd_u64_x(pg_all, z2, add2);
                z3 = svadd_u64_x(pg_all, z3, add3);

                z0 = ReduceModFactor2To1SVE(pg_all, z0, vmod);
                z1 = ReduceModFactor2To1SVE(pg_all, z1, vmod);
                z2 = ReduceModFactor2To1SVE(pg_all, z2, vmod);
                z3 = ReduceModFactor2To1SVE(pg_all, z3, vmod);
            }

            svst1_u64(pg_all, res + i, z0);
            svst1_u64(pg_all, res + i + lanes, z1);
            svst1_u64(pg_all, res + i + 2 * lanes, z2);
            svst1_u64(pg_all, res + i + 3 * lanes, z3);
        }
    
        // Handle the remaining elements with predication.
        while (i < n) {

            const svbool_t pg = svwhilelt_b64(i, n);

            svuint64_t x = svld1_u64(pg, op1 + i);

            // Reduce input to [0, q).
            if constexpr (ModFactor != 1) {
                x = ReduceInputSVE<ModFactor>(pg, x, vmod, v2mod, v4mod);
            }

            // Full 64x64 -> 128-bit product.
            svuint64_t prod_hi;
            svuint64_t prod_lo;

            MulU64ToU128SVE(pg, x, yop, &prod_hi, &prod_lo);

            const svuint64_t c1 =
                ShiftRight128LowPart<Shift>(pg, prod_hi, prod_lo);

            // Barrett quotient approximation.
            const svuint64_t q_hat = svmulh_u64_x(pg, c1, vbarr);

            // Barrett residual: z = prod_lo - q_hat * mod.
            const svuint64_t q_mul = svmul_u64_x(pg, q_hat, vmod);
            svuint64_t z = svsub_u64_x(pg, prod_lo, q_mul);

            // Final correction: [0, 4q) -> [0, q).
            z = ReduceInputSVE<4>(pg, z, vmod, v2mod, v4mod);

            // Add op3 when provided.
            if (op3 != nullptr) {

                svuint64_t add_val = svld1_u64(pg, op3 + i);

                // Reduce addend to [0, q).
                if constexpr (ModFactor != 1) {
                    add_val = ReduceInputSVE<ModFactor>(
                        pg, add_val, vmod, v2mod, v4mod);
                }

                // Modular addition: [0, q) + [0, q) -> [0, 2q) -> [0, q).
                z = svadd_u64_x(pg, z, add_val);
                z = ReduceModFactor2To1SVE(pg, z, vmod);
            }

            svst1_u64(pg, res + i, z);

            i += lanes;
        }
    }


    // Dispatches SVE modular multiply-add to the kernel matching the runtime shift value.
    template <int ModFactor>
    void DispatchEltwiseFMAModSVEKernel(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t barr_factor, uint64_t prod_right_shift) {
    #define AXHEL_DISPATCH_SHIFT(SHIFT_VALUE)                                      \
        case SHIFT_VALUE:                                                          \
            EltwiseFMAModSVEKernel<ModFactor, SHIFT_VALUE>(                        \
                res, op1, op2, op3, n, mod, barr_factor);                          \
            break

        switch(prod_right_shift) {
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
            default: ;
        }
    #undef AXHEL_DISPATCH_SHIFT
    }


    // Performs SVE element-wise modular multiply-add.
    template <int ModFactor>
    void EltwiseFMAModSVE(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod) {
        
        // Barrett params
        constexpr int64_t beta = -2;
        constexpr int64_t alpha = 62;   // alpha - beta = 64

        const uint64_t ceil_log_mod = BitWidth(mod);    // n

        const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

        // Barrett factor mu
        const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

        // Dispatch to the kernel matching the runtime shift value.
        DispatchEltwiseFMAModSVEKernel<ModFactor>(res, op1, op2, op3, n, mod, barr_factor, prod_right_shift);
    }


    // Explicit template instantiations for the supported modulus factors.
    template void EltwiseFMAModSVE<1>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    template void EltwiseFMAModSVE<2>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    template void EltwiseFMAModSVE<4>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    template void EltwiseFMAModSVE<8>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);

}
}

#endif