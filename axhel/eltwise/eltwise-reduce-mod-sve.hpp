// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    void EltwiseReduceModSVE(uint64_t* res, const uint64_t* op, uint64_t n, uint64_t mod, uint64_t in_mod_factor, uint64_t out_mod_factor);

}
}
