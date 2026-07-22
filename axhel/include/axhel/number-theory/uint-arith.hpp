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

    /// @brief Return the Log2 of the number
    inline uint64_t Log2(uint64_t x) noexcept {
        return 63ULL - static_cast<uint64_t>(__builtin_clzll(x));
    }


    /// @brief Return the bit witdth of the number
    inline uint64_t BitWidth(uint64_t x) noexcept {
        return 64ULL - static_cast<uint64_t>(__builtin_clzll(x));
    }


    /// @brief exact scalar 64x64 -> 128 multiplication
    inline void MultiplyUInt64(uint64_t x, uint64_t y, uint64_t* hi, uint64_t* lo) noexcept {
        const uint128_t prod = static_cast<uint128_t>(x) * static_cast<uint128_t>(y);

        *lo = static_cast<uint64_t>(prod);
        *hi = static_cast<uint64_t>(prod >> 64);
    }


    /// @brief High 64 bits of a scalar 64x64 unsigned multiplication.
    inline uint64_t MulHighU64(uint64_t x, uint64_t y) noexcept {
        const uint128_t product = static_cast<uint128_t>(x) * static_cast<uint128_t>(y);

        return static_cast<uint64_t>(product >> 64);
    }


    /// @brief Multiply by a precomputed modular operand floor(multiplier * 2^64 / modulus). The returned result is a lazy modular product.
    inline uint64_t MultiplyUIntModLazy(uint64_t value, uint64_t multiplier, uint64_t multiplier_quotient, uint64_t modulus) noexcept {
        const uint64_t quotient = MulHighU64(value, multiplier_quotient);

        return value * multiplier - quotient * modulus;
    }


    /// @brief Computes the precomputed Shoup quotient  ⌊operand · 2^64 / modulus⌋ for modular multiplication.
    inline uint64_t ComputeShoupQuotient(uint64_t operand, uint64_t modulus) noexcept {
        const uint128_t numerator = static_cast<uint128_t>(operand) << 64;
        
        return static_cast<uint64_t>(numerator / modulus);
    }


} // namespace axhel
} // namespace unipi