// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

class MultiplyFactor {
public:

    /// @brief Compute and store the Barrett factor as floor((operand << bit_shift) / modulus)
    MultiplyFactor(uint64_t operand, uint64_t bit_shift, uint64_t modulus);

    /// @brief Return the operand
    uint64_t Operand() const noexcept;

    /// @brief Return the Barrett factor
    uint64_t BarrettFactor() const noexcept;

private:
    uint64_t _operand;
    uint64_t _barrett_factor;
};

} 
} 