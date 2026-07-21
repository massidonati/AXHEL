// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/ntt/ntt.hpp"

#include <stddef.h>
#include <stdint.h>

namespace unipi {
namespace axhel {

    void NTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers);

    void InverseNTTNegacyclicHarveyLazyNative(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

} // namespace axhel
} // namespace unipi
