// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stddef.h>
#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Operand and Shoup quotient used for modular multiplication.
    struct NTTMultiplyOperand
    {
        uint64_t operand;
        uint64_t quotient;
    };

    /// @brief Forward lazy negacyclic Harvey NTT.
    void NTTNegacyclicHarveyLazy(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief Forward negacyclic Harvey NTT.
    void NTTNegacyclicHarvey(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief Inverse lazy negacyclic Harvey NTT.
    void InverseNTTNegacyclicHarveyLazy(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

    /// @brief Inverse negacyclic Harvey NTT.
    void InverseNTTNegacyclicHarvey(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

} // namespace axhel
} // namespace unipi
