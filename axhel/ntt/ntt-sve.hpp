// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stddef.h>
#include <stdint.h>
#include "axhel/util/defines.hpp"
#include "axhel/ntt/ntt.hpp"


#ifdef AXHEL_HAS_SVE

namespace unipi {
namespace axhel {

    /// @brief Computes the in-place forward negacyclic Harvey NTT in lazy form using SVE.
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    void NTTNegacyclicHarveyLazySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief Computes the in-place normalized forward negacyclic Harvey NTT using SVE.
    /// @post Each output coefficient is in the range [0, modulus).
    void NTTNegacyclicHarveySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    /// @brief Computes the in-place inverse negacyclic Harvey NTT in lazy form using SVE.
    /// @post Each output coefficient is in the range [0, 2 * modulus).
    void InverseNTTNegacyclicHarveyLazySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

    /// @brief Computes the in-place normalized inverse negacyclic Harvey NTT using SVE.
    /// @post Each output coefficient is in the range [0, modulus).
    void InverseNTTNegacyclicHarveySVE(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

} // namespace axhel
} // namespace unipi

#endif