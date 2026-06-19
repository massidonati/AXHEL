// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Compute fused multiply-add (arg1 * arg2 + arg3) mod mod element-wise
    void EltwiseFMAMod(uint64_t* res, const uint64_t* op1, uint64_t op2, const uint64_t* op3, uint64_t n, uint64_t mod, uint64_t mod_factor);

}
}