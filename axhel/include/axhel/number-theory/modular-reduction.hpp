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



    /// @brief Conditionally subtract a scalar bound.
    /// Computes x >= bound ? x - bound : x.
    inline uint64_t ConditionalSubtractNative(uint64_t x, uint64_t bound) noexcept {
        const uint64_t reduced = x - bound;
        return (x >= bound) ? reduced : x;
    }


    /// @brief Reduce a value from [0, 8q) to [0, 4q).
    inline uint64_t ReduceModFactor8To4Native(uint64_t x, uint64_t four_times_mod) noexcept {
        return ConditionalSubtractNative(x, four_times_mod);
    }


    /// @brief Reduce a value from [0, 4q) to [0, 2q).
    inline uint64_t ReduceModFactor4To2Native(uint64_t x, uint64_t twice_mod) noexcept {
        return ConditionalSubtractNative(x, twice_mod);
    }


    /// @brief Reduce a value from [0, 2q) to [0, q).
    inline uint64_t ReduceModFactor2To1Native(uint64_t x, uint64_t modulus) noexcept {
        return ConditionalSubtractNative(x, modulus);
    }


    /// @brief Reduce input value from [0, ModFactor*q) to [0, q). 
    /// For ModFactor = 1 this is intentionally a no-op.
    template <int ModFactor>
    inline uint64_t ReduceInputNative(uint64_t x, uint64_t mod) noexcept {
        
        // [0, 8q) -> [0, 4q).
        if constexpr (ModFactor == 8) {
            x = ReduceModFactor8To4Native(x, 4 * mod);
        }
        
         // [0, 4q) -> [0, 2q).
        if constexpr (ModFactor >= 4) {
            x = ReduceModFactor4To2Native(x, 2 * mod);
        }

         // [0, 2q) -> [0, q).
        if constexpr (ModFactor >= 2) {
            x = ReduceModFactor2To1Native(x, mod);
        }

        return x;
    }



#ifdef AXHEL_HAS_SVE

    /// @brief Conditionally subtracts a lane-wise bound.
    /// For each active lane, computes x >= bound ? x - bound : x.
    inline svuint64_t ConditionalSubtractSVE(svbool_t pg, svuint64_t x, svuint64_t bound) noexcept {
        const svuint64_t reduced = svsub_u64_x(pg, x, bound);
        return svmin_u64_m(pg, x, reduced);
    }


    /// @brief Reduce lane-wise values from [0, 8q) to [0, 4q).
    inline svuint64_t ReduceModFactor8To4SVE(svbool_t pg, svuint64_t x, svuint64_t four_times_mod) noexcept {
        return ConditionalSubtractSVE(pg, x, four_times_mod);
    }


    /// @brief Reduce lane-wise values from [0, 4q) to [0, 2q).
    inline svuint64_t ReduceModFactor4To2SVE(svbool_t pg, svuint64_t x, svuint64_t twice_mod) noexcept {
        return ConditionalSubtractSVE(pg, x, twice_mod);
    }


    /// @brief Reduce lane-wise values from [0, 2q) to [0, q).
    inline svuint64_t ReduceModFactor2To1SVE(svbool_t pg, svuint64_t x, svuint64_t modulus) noexcept {
        return ConditionalSubtractSVE(pg, x, modulus);
    }


    /// @brief Reduce input value from [0, ModFactor*q) to [0, q).
    /// For ModFactor = 1 this is intentionally a no-op.
    template <int ModFactor>
    inline svuint64_t ReduceInputSVE(svbool_t pg, svuint64_t x, svuint64_t vmod, svuint64_t v2mod, svuint64_t v4mod) noexcept {

        // [0, 8q) -> [0, 4q).
        if constexpr (ModFactor == 8) {
            x = ReduceModFactor8To4SVE(pg, x, v4mod);
        }

        // [0, 4q) -> [0, 2q).
        if constexpr (ModFactor >= 4) {
            x = ReduceModFactor4To2SVE(pg, x, v2mod);
        }

        // [0, 2q) -> [0, q).
        if constexpr (ModFactor >= 2) {
            x = ReduceModFactor2To1SVE(pg, x, vmod);
        }

        return x;
    }

#endif


}
}