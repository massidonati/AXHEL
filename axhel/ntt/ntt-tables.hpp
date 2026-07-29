// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/ntt/ntt.hpp"
#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace unipi {
namespace axhel {
namespace detail {

    struct NTTTableData {
        std::vector<NTTMultiplyOperand> root_powers;
        std::vector<NTTMultiplyOperand> inv_root_powers;
        NTTMultiplyOperand inv_degree_modulo;
    };

    NTTTableData GenerateNTTTables(std::size_t degree, std::size_t coeff_count_power, uint64_t modulus, uint64_t root_of_unity);

} // namespace detail
} // namespace axhel
} // namespace unipi