// Copyright (C) 2020 Intel Corporation
// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

#include "axhel/util/defines.hpp"

#ifdef AXHEL_HAS_SVE
#include <arm_sve.h>
#endif

namespace unipi {
namespace axhel {

    /// @brief Reduce a value from [0, 4q) to [0, 2q).
    inline uint64_t ReduceModFactor4To2Native( uint64_t x, uint64_t twice_mod) noexcept {
        if (x >= twice_mod) {
            x -= twice_mod;
        }
        return x;
    }


    /// @brief Reduce a value from [0, 2q) to [0, q).
    inline uint64_t ReduceModFactor2To1Native(uint64_t x, uint64_t modulus) noexcept {
        if (x >= modulus) {
            x -= modulus;
        }
        return x;
    }


    /// @brief Reduce input value from [0, ModFactor*q) to [0, q). For ModFactor = 1 this is intentionally a no-op.
    template <int ModFactor>
    inline uint64_t ReduceInputNative(uint64_t x, uint64_t mod) noexcept {
        
        // if x can be in [0, 4q), subtract 2q
        if constexpr (ModFactor == 4) {
            x = ReduceModFactor4To2Native(x, 2 * mod);
        }

         // if x can be in [0, 2q), subtract q
        if constexpr (ModFactor >= 2) {
        x = ReduceModFactor2To1Native(x, mod);
        }

        return x;
    }

    #ifdef AXHEL_HAS_SVE


    /// @brief Reduce lane-wise values from [0, 4q) to [0, 2q).
    inline svuint64_t ReduceModFactor4To2SVE(svbool_t pg, svuint64_t x, svuint64_t twice_mod) noexcept {
        const svbool_t ge_twice_mod = svcmpge_u64(pg, x, twice_mod);
        return svsub_u64_m(ge_twice_mod, x, twice_mod);
    }

    /// @brief Reduce lane-wise values from [0, 2q) to [0, q).
    inline svuint64_t ReduceModFactor2To1SVE(svbool_t pg, svuint64_t x, svuint64_t modulus) noexcept {
        const svbool_t ge_modulus = svcmpge_u64(pg, x, modulus);
        return svsub_u64_m(ge_modulus, x, modulus);
    }


    /// @brief Reduce input value from [0, ModFactor*q) to [0, q). For ModFactor = 1 this is intentionally a no-op. 
    template <int ModFactor>
    inline svuint64_t ReduceInputSVE(svbool_t pg, svuint64_t x, svuint64_t vmod, svuint64_t v2mod) {

        // if x can be in [0, 4q), subtract 2q
        if constexpr (ModFactor == 4) {
            x = ReduceModFactor4To2SVE(pg, x, v2mod);
        }

         // if x can be in [0, 2q), subtract q
        if constexpr (ModFactor >= 2) {
            x = ReduceModFactor2To1SVE(pg, x, vmod);
        }

        return x;
    }
    #endif

}
}