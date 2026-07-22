// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/ntt/ntt.hpp"

#include <stddef.h>
#include <stdint.h>

#ifdef AXHEL_HAS_SVE

namespace unipi {
namespace axhel {

    /// @brief In-place forward negacyclic Harvey NTT in lazy form using SVE.
    void NTTNegacyclicHarveyLazySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief In-place inverse negacyclic Harvey NTT in lazy form using SVE.
    void InverseNTTNegacyclicHarveyLazySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

} // namespace axhel
} // namespace unipi

#endif