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


    /// @brief Native forward Harvey butterfly. 
    ///
    /// Input:
    ///
    ///   x in [0, 4q)
    ///   y in [0, 4q)
    ///
    /// Output:
    ///
    ///   x in [0, 4q)
    ///   y in [0, 4q)
    ///
    inline void ForwardButterflyNative(uint64_t &x, uint64_t &y, const NTTMultiplyOperand &root, uint64_t modulus, uint64_t twice_modulus) noexcept {
        /*
        * Bring x from [0, 4q) to [0, 2q)
        */
        const uint64_t guarded_x = ReduceModFactor4To2Native(x, twice_modulus);

        /*
        * Lazy Shoup multiplication: transformed_y = y * root mod q represented in [0, 2q).
        */
        const uint64_t transformed_y = MultiplyUIntModLazy(y, root.operand, root.quotient, modulus);

        /*
        * Both results are in [0, 4q).
        */
        x = guarded_x + transformed_y;

        y = guarded_x + twice_modulus - transformed_y;
    }

    /// @brief Native inverse Harvey butterfly. 
    ///
    /// Input:
    ///
    ///   x in [0, 4q)
    ///   y in [0, 2q)
    ///
    /// Output:
    ///
    ///   x in [0, 4q)
    ///   y in [0, 2q)
    ///
    inline void InverseButterflyNative(uint64_t &x, uint64_t &y, const NTTMultiplyOperand &inv_root, uint64_t modulus, uint64_t twice_modulus) noexcept {
        /*
        * Bring x from [0, 4q) to [0, 2q).
        */
        const uint64_t guarded_x = ReduceModFactor4To2Native(x, twice_modulus);

        /*
        * Sum remains in [0, 4q).
        */
        const uint64_t sum = guarded_x + y;

        /*
        * Difference remains non-negative and is in [0, 4q).
        */
        const uint64_t difference = guarded_x + twice_modulus - y;

        x = sum;

        /*
        * Lazy modular product is in [0, 2q).
        */
        y = MultiplyUIntModLazy(difference, inv_root.operand, inv_root.quotient, modulus);
    }


} // namespace detail
} // namespace axhel
} // namespace unipi