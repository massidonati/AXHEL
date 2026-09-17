// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/eltwise/eltwise-mul-mod.hpp" 
#include "eltwise/eltwise-mul-mod-sve.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "eltwise/eltwise-mul-mod-native.hpp"
#include "axhel/util/compiler.hpp"
#include "axhel/util/debug.hpp"

namespace unipi {
namespace axhel {

    // Multiplies two vectors elementwise with modular reduction.
    template <int ModFactor>
    void EltwiseMulModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {

        // Barrett params.
        constexpr int64_t beta = -2;
        constexpr int64_t alpha = 62; // alpha - beta = 64

        const uint64_t ceil_log_mod = BitWidth(mod);    // n
       
        const uint64_t prod_right_shift = static_cast<uint64_t>(static_cast<int64_t>(ceil_log_mod) + beta);

        // Barrett factor mu.
        const uint64_t barr_factor = MultiplyFactor(uint64_t(1) << (ceil_log_mod + alpha - 64), 64, mod).BarrettFactor();

        AXHEL_UNROLL(4)
        // Scalar loop.
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t prod_hi;
            uint64_t prod_lo;
            uint64_t c2_hi;
            uint64_t c2_lo;

            // Reduce inputs to [0, q).
            const uint64_t x = ReduceInputNative<ModFactor>(op1[i], mod);
            const uint64_t y = ReduceInputNative<ModFactor>(op2[i], mod);

            // Full 64x64 -> 128 product
            MultiplyUInt64(x, y, &prod_hi, &prod_lo);

            // c1 = floor((prod_hi:prod_lo) / 2^(n + beta))
            const uint64_t c1 = (prod_lo >> prod_right_shift) | (prod_hi << (64 - prod_right_shift));

            // c2 = floor(U / 2^{n + beta}) * mu
            MultiplyUInt64(c1, barr_factor, &c2_hi, &c2_lo);

            // q_hat = high64(c1 * barr_factor)
            const uint64_t q_hat = c2_hi;

            // Barrett residual.
            const uint64_t z = prod_lo - q_hat * mod;

            // Final correction: [0, 4q) -> [0, q).
            res[i] = ReduceInputNative<4>(z, mod);

        }
    }

    
    // Dispatches element-wise modular multiplication to the available implementation.
    template <int ModFactor>
    inline void EltwiseMulModDispatch(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod) {
        #ifdef AXHEL_HAS_SVE
        AXHEL_LOG("EltwiseMulMod -> SVE" << ", mod_factor=" << ModFactor << ", n=" << n);
        EltwiseMulModSVE<ModFactor>(res, op1, op2, n, mod);
        #else
        AXHEL_LOG("EltwiseMulMod -> native" << ", mod_factor=" << ModFactor << ", n=" << n);
        EltwiseMulModNative<ModFactor>(res, op1, op2, n, mod);
        #endif
    }


    // Dispatches element-wise modular multiplication according to the input modulus factor. 
    // mod_factor must be 1, 2, or 4.
    void EltwiseMulMod(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod, uint64_t mod_factor) {

        switch (mod_factor) {
            case 1:
                EltwiseMulModDispatch<1>(res, op1, op2, n, mod);
                break;
            case 2:
                EltwiseMulModDispatch<2>(res, op1, op2, n, mod);
                break;
            case 4:
                EltwiseMulModDispatch<4>(res, op1, op2, n, mod);
                break;
        }
    }

    // Explicit template instantiations for the supported modulus factors.
    template void EltwiseMulModNative<1>(uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);
    template void EltwiseMulModNative<2>( uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);
    template void EltwiseMulModNative<4>( uint64_t*, const uint64_t*, const uint64_t*, uint64_t, uint64_t);

}
}