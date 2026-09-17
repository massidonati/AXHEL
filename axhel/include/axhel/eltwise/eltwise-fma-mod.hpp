// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Computes element-wise modular multiply-add: (op1 * op2 + op3) mod mod.
    /// @pre Each element of op1 and the scalar op2 must be in the range [0, mod_factor * mod).
    /// @pre If op3 is not null, each element of op3 must be in the range [0, mod_factor * mod).
    /// @pre mod_factor must be 1, 2, 4, or 8.
    void EltwiseFMAMod(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t mod_factor);

}
}