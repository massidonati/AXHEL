// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    /// @brief Subtracts two vectors elementwise with modular reduction.
    /// @pre Each element of op1 and op2 must be in the range [0, mod).
    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod);

    /// @brief Subtracts a vector and scalar elementwise with modular reduction.
    /// @pre Each element of op1 and the scalar op2 must be in the range [0, mod).
    void EltwiseSubModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod);

}
}