// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "ntt/ntt-native.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "ntt/ntt-internal.hpp"
#include "axhel/util/compiler.hpp"
#include "axhel/util/debug.hpp"
#include <stddef.h>
#include <stdint.h>

namespace unipi {
namespace axhel {

    /*
    * Native forward lazy NTT.
    */
    void NTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers) {
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus =  modulus << 1;

        /*
        * Cooley-Tukey stage organization.
        *
        * At every stage:
        *
        *   m   = number of blocks
        *   gap = distance between butterfly operands
        */
        size_t m = 1;
        size_t gap = coeff_count >> 1;

        while (m < coeff_count) {
            for (size_t i = 0; i < m; ++i) {
                /*
                * SEAL's forward root table is stored in bit-reversed
                * order. At stage m, roots are at indices:
                *
                *   m, m+1, ..., 2m-1
                */
                const NTTMultiplyOperand &root = root_powers[m + i];
                const size_t block_start = i * (gap << 1);
                uint64_t *x_ptr = operand + block_start;
                uint64_t *y_ptr = x_ptr + gap;

                AXHEL_UNROLL(4)
                for (size_t j = 0; j < gap; ++j) {
                    detail::ForwardButterflyNative(x_ptr[j], y_ptr[j], root, modulus, twice_modulus);
                }
            }

            m <<= 1;
            gap >>= 1;
        }
    }

    /*
    * Native inverse lazy NTT.
    */
    void InverseNTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo) {
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;

        /*
        * Gentleman-Sande inverse stages.
        *
        * The inverse root table is already stored by SEAL in the exact
        * traversal order required here. Therefore roots must be consumed
        * sequentially.
        */
        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        /*
        * Handle every stage except the final stage.
        *
        * The final stage is specialized because it also incorporates
        * multiplication by n^{-1}.
        */
        while (m > 1) {
            for (size_t i = 0; i < m; ++i) {
                const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];
                const size_t block_start = i * (gap << 1);
                uint64_t *x_ptr = operand + block_start;
                uint64_t *y_ptr =  x_ptr + gap;

                AXHEL_UNROLL(4)
                for (size_t j = 0; j < gap; ++j) {
                    detail::InverseButterflyNative(x_ptr[j], y_ptr[j], inv_root, modulus, twice_modulus);
                }
            }

            m >>= 1;
            gap <<= 1;
        }

        /*
        * At this point root_index identifies the root of the final
        * inverse stage.
        *
        * For N coefficients, this is normally index N-1.
        */
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        /*
        * Construct the combined multiplier:
        *
        *   final_inv_root * n^{-1} mod q
        *
        * SEAL performs the same combination through mul_root_scalar.
        */
        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{
            scaled_root_operand,
            ComputeShoupQuotient(scaled_root_operand, modulus)
        };

        /*
        * Final inverse stage.
        *
        * At this point:
        *
        *   m   = 1
        *   gap = coeff_count / 2
        *
        * Every upper result is multiplied by n^{-1}.
        * Every lower result is multiplied by final_root * n^{-1}.
        */
        uint64_t *x_ptr = operand;

        uint64_t *y_ptr = operand + gap;

        AXHEL_UNROLL(4)
        for (size_t j = 0; j < gap; ++j) {
            const uint64_t guarded_x = ReduceModFactor4To2Native(x_ptr[j], twice_modulus);
            const uint64_t y = y_ptr[j];
            const uint64_t sum = guarded_x + y;
            const uint64_t difference = guarded_x + twice_modulus - y;

            /*
            * Both lazy products are returned in [0, 2q).
            */
            x_ptr[j] = MultiplyUIntModLazy(sum, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

            y_ptr[j] = MultiplyUIntModLazy( difference, scaled_root.operand, scaled_root.quotient, modulus);
        }
    }


} // namespace axhel
} // namespace unipi    