// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-reduce-mod.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "eltwise/eltwise-reduce-mod-native.hpp"
#include "eltwise/eltwise-reduce-mod-sve.hpp"
#include "axhel/util/compiler.hpp"

namespace unipi {
namespace axhel {

    inline uint64_t BarrettReduceUInt64Native(uint64_t op, uint64_t mod, uint64_t barr_factor, uint64_t prod_right_shift) {
        uint64_t c2_hi;
        uint64_t c2_lo;

        const uint64_t c1 = op >> prod_right_shift;
        MultiplyUInt64(c1, barr_factor, &c2_hi, &c2_lo);
        const uint64_t q_hat = c2_hi;
        const uint64_t z = op - q_hat * mod;

        return ReduceInputNative<4>(z, mod);
    }


    void EltwiseReduceModNative(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor) {

        if(in_mod_factor == 1) {
            for(uint64_t i=0; i<n; ++i){
                *res = *op;
                ++op;
                ++res;
            }
            return;
        }

        if(out_mod_factor == 2 && in_mod_factor <= 2) {
            for(uint64_t i=0; i<n; ++i){
                *res = *op;
                ++op;
                ++res;
            }
            return;
        }

        if(in_mod_factor == 2 || in_mod_factor == 4) {
            const uint64_t target = out_mod_factor * mod;

            for(uint64_t i=0; i<n; ++i){
                uint64_t res_val = *op;

                while(res_val >= target) {
                    res_val -= mod;
                }

                *res = res_val;

                ++op;
                ++res;
            }

            return;
        }

        if(in_mod_factor == mod) {
            constexpr int64_t beta = -2;
            constexpr int64_t alpha = 62;

            const uint64_t ceil_log_mod = BitWidth(mod);
            const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

            const uint64_t barr_factor =  MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

            AXHEL_UNROLL(4)
            for(uint64_t i=0; i<n; ++i){
                *res = BarrettReduceUInt64Native(*op, mod, barr_factor, prod_right_shift);

                ++op;
                ++res;
            }

            return;
        }
    }


    void EltwiseReduceMod(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor) {
        #ifdef AXHEL_HAS_SVE
        EltwiseReduceModSVE(res,op,n,mod,in_mod_factor,out_mod_factor);
        #else
        EltwiseReduceModNative(res,op,n,mod,in_mod_factor,out_mod_factor);
        #endif
    }

}
}