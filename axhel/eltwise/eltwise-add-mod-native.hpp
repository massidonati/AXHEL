// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace unipi {
namespace axhel {

    /// @brief Adds two vectors elementwise with modular reduction
    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, const uint64_t* op2, uint64_t n, uint64_t mod);

    /// @brief Adds a vector and scalar elementwise with modular reduction
    void EltwiseAddModNative(uint64_t* res, const uint64_t* op1, uint64_t op2, uint64_t n, uint64_t mod);

}
}