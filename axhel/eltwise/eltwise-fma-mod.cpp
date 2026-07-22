// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-fma-mod.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "eltwise/eltwise-fma-mod-native.hpp"

#ifdef AXHEL_HAS_SVE
#include "eltwise/eltwise-fma-mod-sve.hpp"
#endif

namespace unipi {
namespace axhel {

    template <int ModFactor>
    void EltwiseFMAModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod) {
       
        //TODO: check

         // Barrett params
        constexpr int64_t beta = -2;
        constexpr int64_t alpha = 62; // alpha - beta = 64

        const uint64_t ceil_log_mod = BitWidth(mod); // n

        const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

        // Barrett factor mu
        const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

        // modular reduction of input
        const uint64_t op2_red = ReduceInputNative<ModFactor>(op2, mod);

        for(uint64_t i=0; i<n; ++i) {
            uint64_t prod_hi;
            uint64_t prod_lo;
            uint64_t c2_hi;
            uint64_t c2_lo;

            // modular reduction of input
            const uint64_t x = ReduceInputNative<ModFactor>(op1[i], mod);

            // full 64x64 -> 128 product
            MultiplyUInt64(x, op2_red, &prod_hi, &prod_lo);

             // c1 = floor((prod_hi:prod_lo) / 2^(n + beta))
            const uint64_t c1 = (prod_lo >> prod_right_shift) | (prod_hi << (64 - prod_right_shift));

            // c2 = floor(U / 2^{n + beta}) * mu
            MultiplyUInt64(c1, barr_factor, &c2_hi, &c2_lo);

            // q_hat = high64(c1 * barr_factor)
            const uint64_t q_hat = c2_hi;

            // only the low 64 bits are required here
            const uint64_t z = prod_lo - q_hat * mod;

            // final correction to [0, q)
            uint64_t res_val = ReduceInputNative<4>(z, mod);

            // conditional addition
            if(op3 != nullptr) {
                uint64_t add_val = ReduceInputNative<ModFactor>(op3[i], mod);

                res_val += add_val;

                if(res_val >= mod) {
                    res_val -= mod;
                }
            }

            res[i] = res_val;
        }
    }


    template <int ModFactor>
    inline void EltwiseFMAModDispatch(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod) {
        #ifdef AXHEL_HAS_SVE
            EltwiseFMAModSVE<ModFactor>(res, op1, op2, op3, n, mod);
        #else
            EltwiseFMAModNative<ModFactor>(res, op1, op2, op3, n, mod);
        #endif
    }


    void EltwiseFMAMod(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t mod_factor) {
        
        //TODO: check

        switch(mod_factor) {
            case 1:
                EltwiseFMAModDispatch<1>(res, op1, op2, op3, n, mod);
                break;
            case 2:
                EltwiseFMAModDispatch<2>(res, op1, op2, op3, n, mod);
                break;
            case 4:
                EltwiseFMAModDispatch<4>(res, op1, op2, op3, n, mod);
                break;
        }
    }


    template void EltwiseFMAModNative<1>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    
    template void EltwiseFMAModNative<2>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    
    template void EltwiseFMAModNative<4>(uint64_t*, const uint64_t*, uint64_t, const uint64_t*, uint64_t, uint64_t);
    
}
}