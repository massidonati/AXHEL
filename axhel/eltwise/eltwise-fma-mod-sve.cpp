// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>

#include "eltwise/eltwise-fma-mod-sve.hpp"

#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/sve-arith.hpp"
#include "axhel/number-theory/uint-arith.hpp"

#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>

namespace unipi {
namespace axhel {

    //@brief SVE specialized internal kernel 
    template <int ModFactor, uint64_t Shift>
    void EltwiseFMAModSVEKernel(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t barr_factor) {
        const svuint64_t vmod = svdup_n_u64(mod);
        const svuint64_t v2mod = svdup_n_u64(2 * mod);
        const svuint64_t vbarr = svdup_n_u64(barr_factor);

        uint64_t op2_red = ReduceInputNative<ModFactor>(op2, mod);
        const svuint64_t yop = svdup_n_u64(op2_red);

        const uint64_t lanes = svcntd();

        #pragma unroll 4
        for(uint64_t i=0; i<n; i+=lanes) {
            svbool_t pg = svwhilelt_b64(i, n);

            svuint64_t x = svld1_u64(pg, op1 + i);

            if constexpr (ModFactor != 1) {
                x = ReduceInputSVE<ModFactor>(pg, x, vmod, v2mod);
            }

            svuint64_t prod_hi;
            svuint64_t prod_lo;
            MulU64ToU128SVE(pg, x, yop, &prod_hi, &prod_lo);

            svuint64_t c1 = ShiftRight128LowPart<Shift>(pg, prod_hi, prod_lo);
            svuint64_t q_hat = svmulh_u64_x(pg, c1, vbarr);

            svuint64_t q_mul = svmul_u64_x(pg, q_hat, vmod);
            svuint64_t z = svsub_u64_x(pg, prod_lo, q_mul);

            z = ReduceInputSVE<4>(pg, z, vmod, v2mod);

            if(op3 != nullptr) {
                svuint64_t add_val = svld1_u64(pg, op3 + i);

                if constexpr (ModFactor != 1) {
                    add_val = ReduceInputSVE<ModFactor>(pg, add_val, vmod, v2mod);
                }

                z = svadd_u64_x(pg, z, add_val);

                svbool_t ge_mod = svcmpge_u64(pg, z, vmod);
                z = svsub_u64_m(ge_mod, z, vmod);
            }

            svst1_u64(pg, res + i, z);
        }
    }


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


    template <int ModFactor>
    void EltwiseFMAModSVE(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod) {
        
        //TODO: check
        
        // Barrett params
        constexpr int64_t beta = -2;
        constexpr int64_t alpha = 62;   // alpha - beta = 64

        const uint64_t ceil_log_mod = BitWidth(mod);    // n

        const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

        // Barrett factor mu
        const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

        // specialized dispatch 
        DispatchEltwiseFMAModSVEKernel<ModFactor>(res, op1, op2, op3, n, mod, barr_factor, prod_right_shift);
    }


    template void EltwiseFMAModSVE<1>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    
    template void EltwiseFMAModSVE<2>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    
    template void EltwiseFMAModSVE<4>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);

}
}

#endif