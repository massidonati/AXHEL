// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "eltwise/eltwise-mul-mod-sve2.hpp"

#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/sve-arith.hpp"

#include <arm_sve.h>
#include <stdint.h>

#ifdef AXHEL_HAS_SVE2

namespace unipi {
namespace axhel {


    //@brief SVE2 specialized internal kernel 
    template <int ModFactor, uint64_t Shift>
    void EltwiseMulModSVE2Kernel(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod, uint64_t barr_factor) {
        const svuint64_t vmod = svdup_n_u64(mod);
        const svuint64_t v2mod = svdup_n_u64(2 * mod);
        const svuint64_t vbarr = svdup_n_u64(barr_factor);

        const uint64_t lanes = svcntd();

        #pragma unroll 4
        for (uint64_t i = 0; i < n; i += lanes) {
            svbool_t pg = svwhilelt_b64(i, n);

            svuint64_t x = svld1_u64(pg, op1 + i);
            svuint64_t y = svld1_u64(pg, op2 + i);

             // normalize lazy inputs only when required by ModFactor
            if constexpr (ModFactor != 1) {
                x = ReduceInputSVE<ModFactor>(pg, x, vmod, v2mod);
                y = ReduceInputSVE<ModFactor>(pg, y, vmod, v2mod);
            }

            // Multiply 64x64 -> 128 bit, hi and lo parts
            svuint64_t prod_lo = svmul_u64_x(pg, x, y);
            svuint64_t prod_hi = svmulh_u64_x(pg, x, y);

            // Barrett: c1 = (prod_hi:prod_lo) >> Shift
            svuint64_t c1 = ShiftRight128LowPart<Shift>(pg, prod_hi, prod_lo);

            // q_hat = high64(c1 * barr_factor)
            svuint64_t q_hat = svmulh_u64_x(pg, c1, vbarr);

            // z = prod_lo - q_hat * modulus
            svuint64_t q_mul = svmul_u64_x(pg, q_hat, vmod);
            svuint64_t z = svsub_u64_x(pg, prod_lo, q_mul);

            // final correction to [0, q)
            z = ReduceInputSVE<4>(pg, z, vmod, v2mod);

            //TODO: valutare se scrivere istruzioni inline

            svst1_u64(pg, res + i, z);
        }
    }


    // @brief out-of-loop dispatcher to select template version with Shift known at compile-time 
    template <int ModFactor>
    void DispatchEltwiseMulModSVE2Kernel(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod, uint64_t barr_factor, uint64_t prod_right_shift) {
    #define AXHEL_DISPATCH_SHIFT(SHIFT_VALUE)                                      \
        case SHIFT_VALUE:                                                          \
            EltwiseMulModSVE2Kernel<ModFactor, SHIFT_VALUE>(                       \
                res, op1, op2, n, mod, barr_factor);                               \
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
        default: ; //TODO error
        }
    #undef AXHEL_DISPATCH_SHIFT
    }


    // @brief SVE2-optimized multiply
    template <int ModFactor>
    void EltwiseMulModSVE2(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        
        //TODO: check

        // Barrett params
        constexpr int64_t beta = -2;
        constexpr int64_t alpha = 62;   // alpha - beta = 64

        //const uint64_t ceil_log_mod = Log2(modulus) + 1;  // n
        const uint64_t ceil_log_mod = BitWidth(mod);    // n
    
        const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

        // Barrett factor mu
        const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();  

        // specialized dispatch 
        DispatchEltwiseMulModSVE2Kernel<ModFactor>(res, op1, op2, n, mod, barr_factor, prod_right_shift);

    }

    template void EltwiseMulModSVE2<1>(uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);

    template void EltwiseMulModSVE2<2>(uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);

    template void EltwiseMulModSVE2<4>(uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);

}
}

#endif
