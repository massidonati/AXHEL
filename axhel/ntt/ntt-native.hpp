// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/ntt/ntt.hpp"

#include <stddef.h>
#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Computes the in-place forward negacyclic Harvey NTT in lazy form.
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    void NTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief Computes the in-place inverse negacyclic Harvey NTT in lazy form.
    /// @post Each output coefficient is in the range [0, 2 * modulus).
    void InverseNTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

} // namespace axhel
} // namespace unipi
