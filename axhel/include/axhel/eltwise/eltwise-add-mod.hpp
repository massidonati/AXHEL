// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Adds two vectors elementwise with modular reduction.
    /// @pre Each element of operand1 and operand2 must be in the range [0, modulus).
    void EltwiseAddMod(uint64_t* result, const uint64_t* operand1, const uint64_t* operand2, uint64_t n, uint64_t modulus);

    /// @brief Adds a vector and scalar elementwise with modular reduction.
    /// @pre Each element of operand1 and the scalar operand2 must be in the range [0, modulus).
    void EltwiseAddMod(uint64_t* result, const uint64_t* operand1, uint64_t operand2, uint64_t n, uint64_t modulus);

} 
} 