// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include "axhel/util/defines.hpp"

#ifdef AXHEL_HAS_SVE
#include <arm_sve.h>

namespace unipi {
namespace axhel {
namespace test_hooks {

// These functions are defined in axhel/ntt/ntt-sve.cpp only when
// AXHEL_TESTING is enabled. They expose the real anonymous-namespace
// butterfly implementations without making them part of the installed API.
void ForwardButterflySVE(
    svbool_t pg,
    svuint64_t x,
    svuint64_t y,
    svuint64_t root,
    svuint64_t root_quotient,
    svuint64_t modulus,
    svuint64_t twice_modulus,
    svuint64_t& result_x,
    svuint64_t& result_y) noexcept;

void InverseButterflySVE(
    svbool_t pg,
    svuint64_t x,
    svuint64_t y,
    svuint64_t inv_root,
    svuint64_t inv_root_quotient,
    svuint64_t modulus,
    svuint64_t twice_modulus,
    svuint64_t& result_x,
    svuint64_t& result_y) noexcept;

} // namespace test_hooks
} // namespace axhel
} // namespace unipi
#endif
