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

    /// @brief Computes the native in-place forward negacyclic Harvey NTT in lazy form.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
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

        // Process all forward NTT stages.
        while (m < coeff_count) {
            for (size_t i = 0; i < m; ++i) {

                // Forward roots are stored in bit-reversed order.
                const NTTMultiplyOperand &root = root_powers[m + i];

                const size_t block_start = i * (gap << 1);
                uint64_t *x_ptr = operand + block_start;
                uint64_t *y_ptr = x_ptr + gap;

                AXHEL_UNROLL(4)
                // Apply the forward Harvey butterflies for the current block.
                for (size_t j = 0; j < gap; ++j) {
                    detail::ForwardButterflyNative(x_ptr[j], y_ptr[j], root, modulus, twice_modulus);
                }
            }

            m <<= 1;
            gap >>= 1;
        }
    }


    /// @brief Computes the native in-place inverse negacyclic Harvey NTT in lazy form.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).
    void InverseNTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo) {
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        /*
        * Gentleman-Sande stage organization.
        *
        * At every stage:
        *
        *   m   = number of blocks
        *   gap = distance between butterfly operands
        */
        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        /*
         * Process all stages except the final inverse stage.
         * The final stage also incorporates multiplication by n^{-1}.
         */
        while (m > 1) {
            for (size_t i = 0; i < m; ++i) {

                 // Inverse roots are consumed sequentially.
                const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];

                const size_t block_start = i * (gap << 1);
                uint64_t *x_ptr = operand + block_start;
                uint64_t *y_ptr =  x_ptr + gap;

                AXHEL_UNROLL(4)
                // Apply the inverse Harvey butterflies for the current block.
                for (size_t j = 0; j < gap; ++j) {
                    detail::InverseButterflyNative(x_ptr[j], y_ptr[j], inv_root, modulus, twice_modulus);
                }
            }

            m >>= 1;
            gap <<= 1;
        }

        // Root used by the final inverse stage.
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        /*
         * Combine the final inverse root with n^{-1}:
         *
         *   scaled_root = final_inv_root * n^{-1} mod q
         */
        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{ scaled_root_operand, ComputeShoupQuotient(scaled_root_operand, modulus) };

        /*
        * Final inverse stage.
        *
        * At this point:
        *
        *   m   = 1
        *   gap = coeff_count / 2
        *
        * The upper branch is multiplied by n^{-1}.
        * The lower branch is multiplied by scaled_root.
        */
        uint64_t *x_ptr = operand;

        uint64_t *y_ptr = operand + gap;

        AXHEL_UNROLL(4)
        for (size_t j = 0; j < gap; ++j) {
            // Inputs are reduced to the range required by the final stage.
            const uint64_t guarded_x = ReduceModFactor4To2Native(x_ptr[j], twice_modulus);
            const uint64_t y = y_ptr[j];
            const uint64_t sum = guarded_x + y;
            const uint64_t difference = guarded_x + twice_modulus - y;

            // Lazy Shoup products in [0, 2q).
            x_ptr[j] = MultiplyUIntModLazy(sum, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);
            y_ptr[j] = MultiplyUIntModLazy( difference, scaled_root.operand, scaled_root.quotient, modulus);
        }
    }


} // namespace axhel
} // namespace unipi    