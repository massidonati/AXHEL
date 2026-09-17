// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/ntt/ntt.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/uint-arith.hpp"

#include <stdint.h>

namespace unipi {
namespace axhel {
namespace detail {


    /// @brief Performs a native forward Harvey butterfly.
    /// @pre x and y are in the range [0, 4q).
    /// @post x and y are in the range [0, 4q).
    inline void ForwardButterflyNative(uint64_t &x, uint64_t &y, const NTTMultiplyOperand &root, uint64_t modulus, uint64_t twice_modulus) noexcept {
        
        // [0, 4q) -> [0, 2q).
        const uint64_t guarded_x = ReduceModFactor4To2Native(x, twice_modulus);

        // Lazy Shoup product in [0, 2q).
        const uint64_t transformed_y = MultiplyUIntModLazy(y, root.operand, root.quotient, modulus);

        // Results remain in [0, 4q).
        x = guarded_x + transformed_y;
        y = guarded_x + twice_modulus - transformed_y;
    }


    /// @brief Performs a native inverse Harvey butterfly.
    /// @pre x and y are in the range [0, 2q).
    /// @post x and y are in the range [0, 2q).
    inline void InverseButterflyNative(uint64_t &x, uint64_t &y, const NTTMultiplyOperand &inv_root, uint64_t modulus, uint64_t twice_modulus) noexcept {

        const uint64_t u = x;
        const uint64_t v = y;

        // [0, 4q) -> [0, 2q).
        x = ReduceModFactor4To2Native(u + v, twice_modulus);

        const uint64_t difference = u + twice_modulus - v;

        // Lazy Shoup product in [0, 2q).
        y = MultiplyUIntModLazy(difference, inv_root.operand, inv_root.quotient, modulus);
    }


} // namespace detail
} // namespace axhel
} // namespace unipi