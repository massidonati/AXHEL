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
#include "axhel/util/debug.hpp"


namespace unipi {
namespace axhel {

    // Performs Barrett reduction of a 64-bit input using precomputed parameters.
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

        // Input already reduced modulo q.
        if(in_mod_factor == 1) {
            for(uint64_t i=0; i<n; ++i){
                *res = *op;
                ++op;
                ++res;
            }
            return;
        }

        // Input already satisfies the requested [0, 2q) output range.
        if(out_mod_factor == 2 && in_mod_factor <= 2) {
            for(uint64_t i=0; i<n; ++i){
                *res = *op;
                ++op;
                ++res;
            }
            return;
        }

        // Reduce from [0, 2q) to [0, q).
        if(in_mod_factor == 2 && out_mod_factor == 1) {
            for(uint64_t i=0; i<n; ++i) {
                *res = ReduceModFactor2To1Native(*op, mod);

                ++op;
                ++res;
            }

            return;
        }

        // Reduce from [0, 4q) to [0, q).
        if(in_mod_factor == 4 && out_mod_factor == 1) {
            const uint64_t twice_mod = 2 * mod;

            for(uint64_t i=0; i<n; ++i) {
                uint64_t res_val = *op;

                // [0, 4q) -> [0, 2q)
                res_val = ReduceModFactor4To2Native(res_val, twice_mod);

                // [0, 2q) -> [0, q)
                res_val = ReduceModFactor2To1Native(res_val, mod);

                *res = res_val;

                ++op;
                ++res;
            }

            return;
        }

        // Reduce from [0, 4q) to [0, 2q).
        if(in_mod_factor == 4 && out_mod_factor == 2) {
            const uint64_t twice_mod = 2 * mod;

            for(uint64_t i=0; i<n; ++i) {
                uint64_t res_val = *op;

                if(res_val >= twice_mod) {
                    res_val -= twice_mod;
                }

                *res = res_val;

                ++op;
                ++res;
            }

            return;
        }  

        // Barrett reduction for the general input range [0, mod^2).
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


    // Dispatches element-wise modular reduction to the available implementation.
    void EltwiseReduceMod(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor) {
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseReduceMod -> SVE" << ", in_factor=" << in_mod_factor << ", out_factor=" << out_mod_factor << ", n=" << n);
        EltwiseReduceModSVE(res,op,n,mod,in_mod_factor,out_mod_factor);
        #else
        AXHEL_LOG("EltwiseReduceMod -> native" << ", in_factor=" << in_mod_factor << ", out_factor=" << out_mod_factor << ", n=" << n);
        EltwiseReduceModNative(res,op,n,mod,in_mod_factor,out_mod_factor);
        #endif
    }

}
}