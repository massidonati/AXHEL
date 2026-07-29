// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <memory>

namespace unipi {
namespace axhel {

    /// @brief Operand and Shoup quotient used for modular multiplication.
    struct NTTMultiplyOperand
    {
        uint64_t operand;
        uint64_t quotient;
    };


    /// @brief Negacyclic Number Theoretic Transform.
    ///
    /// The object stores or references all precomputed values required by the
    /// forward and inverse transforms.
    class NTT {
    
        public:
        /// @brief Construct an autonomous NTT object.
        ///
        /// AXHEL generates and owns the forward and inverse root tables.
        ///
        /// @param degree Polynomial modulus degree. Must be a power of two.
        /// @param modulus NTT-friendly prime modulus.
        /// @param root_of_unity Primitive 2*degree-th root of unity modulo modulus.
        NTT(uint64_t degree, uint64_t modulus, uint64_t root_of_unity);

        /// @brief Construct an NTT object using externally owned tables.
        ///
        /// No table is copied or generated. The caller must guarantee that the
        /// tables remain alive for the entire lifetime of this NTT object.
        NTT(uint64_t degree, uint64_t modulus, const NTTMultiplyOperand *root_powers, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo);

        ~NTT();

        NTT(const NTT &) noexcept;
        NTT &operator=(const NTT &) noexcept;

        NTT(NTT &&) noexcept;
        NTT &operator=(NTT &&) noexcept;

        /// @brief Compute the forward negacyclic NTT.
        ///
        /// result and operand may point to the same buffer.
        void ComputeForward(uint64_t *result, const uint64_t *operand, uint64_t input_mod_factor = 1, uint64_t output_mod_factor = 1) const;

        /// @brief Compute the inverse negacyclic NTT.
        ///
        /// result and operand may point to the same buffer.
        void ComputeInverse(uint64_t *result, const uint64_t *operand, uint64_t input_mod_factor = 1, uint64_t output_mod_factor = 1) const;

        uint64_t Degree() const noexcept;
        uint64_t Modulus() const noexcept;

    private:
        class Impl;

        std::shared_ptr<Impl> impl_;
    };

} // namespace axhel
} // namespace unipi
