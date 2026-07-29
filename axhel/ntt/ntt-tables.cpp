// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "ntt/ntt-tables.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace unipi {
namespace axhel {
namespace detail {

namespace {

    inline uint64_t MultiplyMod(uint64_t operand1, uint64_t operand2, uint64_t modulus) noexcept {
        const uint128_t product = static_cast<uint128_t>(operand1) * static_cast<uint128_t>(operand2);
        return static_cast<uint64_t>(product % modulus);
    }


    uint64_t PowMod(uint64_t base, uint64_t exponent, uint64_t modulus) noexcept {
        uint64_t result = 1;
        base %= modulus;

        while (exponent != 0) {
            if ((exponent & 1ULL) != 0) {
                result = MultiplyMod(result, base, modulus);
            }

            exponent >>= 1;

            if (exponent != 0) {
                base = MultiplyMod(base, base, modulus);
            }
        }

        return result;
    }


    inline NTTMultiplyOperand MakeNTTMultiplyOperand(uint64_t operand, uint64_t modulus) noexcept {
        return NTTMultiplyOperand{
            operand,
            ComputeShoupQuotient(operand, modulus)
        };
    }


    inline std::size_t NextBitReversedIndex(std::size_t current, std::size_t degree) noexcept {
        std::size_t bit = degree >> 1;

        while ((current & bit) != 0) {
            current ^= bit;
            bit >>= 1;
        }

        current ^= bit;

        return current;
    }


    void GenerateSequentialPowers(std::vector<uint64_t>& powers, uint64_t root, uint64_t modulus) {
        const std::size_t degree = powers.size();

        powers[0] = 1;

        for (std::size_t i = 1; i < degree; ++i) {
            powers[i] = MultiplyMod(powers[i - 1], root, modulus);
        }
    }

} // namespace


    NTTTableData GenerateNTTTables(std::size_t degree, std::size_t coeff_count_power, uint64_t modulus, uint64_t root_of_unity) {
        static_cast<void>(coeff_count_power);

        NTTTableData tables;

        tables.root_powers.resize(degree);
        tables.inv_root_powers.resize(degree);

        /*
        * Temporary storage reused for forward and inverse powers.
        *
        * Only one additional degree-element vector is allocated.
        */
        std::vector<uint64_t> sequential_powers(degree);

        /*

        * Generate:
        *
        *   1, root, root^2, ..., root^(degree-1)
        */
        GenerateSequentialPowers(sequential_powers, root_of_unity, modulus);

        /*
        * Store forward roots in bit-reversed order:
        *
        *   root_powers[i] =
        *       root_of_unity ^ ReverseBits(i)
        */
        std::size_t reversed_index = 0;

        for (std::size_t i = 0; i < degree; ++i) {
            const uint64_t root = sequential_powers[reversed_index];

            tables.root_powers[i] = MakeNTTMultiplyOperand(root, modulus);

            if (i + 1 < degree) {
                reversed_index = NextBitReversedIndex(reversed_index, degree);
            }
        }

        const uint64_t inv_root_of_unity = PowMod(root_of_unity, modulus - 2, modulus);

        GenerateSequentialPowers(sequential_powers, inv_root_of_unity, modulus);

        /*
        * Index zero is unused by the inverse Harvey kernel.
        */
        tables.inv_root_powers[0] = NTTMultiplyOperand{0, 0};

        reversed_index = 0;

        for (std::size_t i = 1; i < degree; ++i) {
            const std::size_t exponent = reversed_index + 1;

            const uint64_t inv_root = sequential_powers[exponent];

            tables.inv_root_powers[i] = MakeNTTMultiplyOperand(inv_root, modulus);

            if (i + 1 < degree) {
                reversed_index = NextBitReversedIndex(reversed_index, degree);
            }
        }

        const uint64_t degree_modulo = static_cast<uint64_t>(degree) % modulus;

        const uint64_t inv_degree = PowMod(degree_modulo, modulus - 2, modulus);

        tables.inv_degree_modulo = MakeNTTMultiplyOperand(inv_degree, modulus);

        return tables;

    }

} // namespace detail
} // namespace axhel
} // namespace unipi


