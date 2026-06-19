// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/number-theory/multiply-factor.hpp"
#include "axhel/number-theory/uint-arith.hpp"

namespace unipi {
namespace axhel {

    MultiplyFactor::MultiplyFactor(
        uint64_t operand,
        uint64_t bit_shift,
        uint64_t modulus)
        : _operand(operand)
    {
        uint128_t numerator = static_cast<uint128_t>(operand) << bit_shift;
        _barrett_factor = static_cast<uint64_t>(numerator / modulus);
    }

    uint64_t MultiplyFactor::Operand() const noexcept
    {
        return _operand;
    }

    uint64_t MultiplyFactor::BarrettFactor() const noexcept
    {
        return _barrett_factor;
    }

}
}