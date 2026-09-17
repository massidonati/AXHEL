// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Performs element-wise modular reduction.
    /// @pre Each element of op must be in the range [0, in_mod_factor * mod).
    /// @post Each element of res is in the range [0, out_mod_factor * mod).
    void EltwiseReduceMod(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor);

} 
} 