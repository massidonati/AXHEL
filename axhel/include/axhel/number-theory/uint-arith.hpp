// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

namespace unipi {
namespace axhel {

    #if defined(__SIZEOF_INT128__)
    using uint128_t = unsigned __int128;
    using int128_t = __int128;
    #else
    #error "AXHEL requires compiler support for __int128"     
    #endif

    //@brief Return the Log2 of the number
    inline uint64_t Log2(uint64_t x) noexcept {
        return 63ULL - static_cast<uint64_t>(__builtin_clzll(x));
    }

    //@brief Return the bit witdth of the number
    inline uint64_t BitWidth(uint64_t x) noexcept {
        return 64ULL - static_cast<uint64_t>(__builtin_clzll(x));
    }

    //@brief exact scalar 64x64 -> 128 multiplication
    inline void MultiplyUInt64(uint64_t x, uint64_t y, uint64_t* hi, uint64_t* lo) {
        const uint128_t prod = static_cast<uint128_t>(x) * static_cast<uint128_t>(y);

        *lo = static_cast<uint64_t>(prod);
        *hi = static_cast<uint64_t>(prod >> 64);
    }


} // namespace axhel
} // namespace unipi