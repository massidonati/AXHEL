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

    template <uint64_t Shift>
    inline svuint64_t BarrettReduceUInt64SVEKernel(svbool_t pg, svuint64_t op_seg, uint64_t mod, uint64_t barr_factor) {
        const svuint64_t vmod = svdup_n_u64(mod);
        const svuint64_t v2mod = svdup_n_u64(2 * mod);
        const svuint64_t v4mod = svdup_n_u64(4 * mod);
        const svuint64_t vbarr = svdup_n_u64(barr_factor);
        const svuint64_t zero = svdup_n_u64(0);

        svuint64_t c1 = ShiftRight128LowPart<Shift>(pg, zero, op_seg);

        svuint64_t q_hat = svmulh_u64_x(pg, c1, vbarr); 

        svuint64_t q_mul = svmul_u64_x(pg, q_hat, vmod);
        svuint64_t z = svsub_u64_x(pg, op_seg, q_mul);

        z = ReduceInputSVE<4>(pg, z, vmod, v2mod, v4mod);

        return z;
    }


    void DispatchBarrettReduceUInt64SVEKernel(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t barr_factor, uint64_t prod_right_shift) {
        const uint64_t lanes = svcntd();

    #define AXHEL_DISPATCH_SHIFT(SHIFT_VALUE)                                      \
        case SHIFT_VALUE:                                                          \
            AXHEL_UNROLL(4)                                                        \
            for(uint64_t i=0; i<n; i+=lanes) {                                     \
                svbool_t pg = svwhilelt_b64(i, n);                                 \
                svuint64_t op_seg = svld1_u64(pg, op + i);                         \
                svuint64_t res_seg = BarrettReduceUInt64SVEKernel<SHIFT_VALUE>(    \
                    pg, op_seg, mod, barr_factor);                                 \
                svst1_u64(pg, res + i, res_seg);                                   \
            }                                                                      \
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
            default: ; //TODO error
        }
    #undef AXHEL_DISPATCH_SHIFT
    }


    void EltwiseReduceModSVE(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor) {

        const uint64_t lanes = svcntd();

        if(in_mod_factor == 1) {
            for(uint64_t i=0; i<n; i+=lanes){
                svbool_t pg = svwhilelt_b64(i, n);
                svuint64_t res_seg = svld1_u64(pg, op + i);
                svst1_u64(pg, res + i, res_seg);
            }
            return;
        }

        if(out_mod_factor == 2 && in_mod_factor <= 2) {
            for(uint64_t i=0; i<n; i+=lanes){
                svbool_t pg = svwhilelt_b64(i, n);
                svuint64_t res_seg = svld1_u64(pg, op + i);
                svst1_u64(pg, res + i, res_seg);
            }
            return;
        }

        if(in_mod_factor == 2 || in_mod_factor == 4) {
            const uint64_t target = out_mod_factor * mod;

            for(uint64_t i=0; i<n; i+=lanes){
                svbool_t pg = svwhilelt_b64(i, n);
                svuint64_t res_seg = svld1_u64(pg, op + i);
                svbool_t ge_cond = svcmpge_n_u64(pg, res_seg, target);

                while(svptest_any(pg, ge_cond)) {
                    res_seg = svsub_n_u64_m(ge_cond, res_seg, mod);
                    ge_cond = svcmpge_n_u64(pg, res_seg, target);
                }

                svst1_u64(pg, res + i, res_seg);
            }

            return;
        }

        if(in_mod_factor == mod) {
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