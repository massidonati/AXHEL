// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "axhel/ntt/ntt.hpp"
#include "axhel/eltwise/eltwise-reduce-mod.hpp"
#include "ntt/ntt-native.hpp"
#include "ntt/ntt-tables.hpp"
#include "axhel/util/debug.hpp"
#include <stddef.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>

#ifdef AXHEL_HAS_SVE
#include "ntt/ntt-sve.hpp"
#endif


namespace unipi {
namespace axhel {


    class NTT::Impl {
    public:
        uint64_t degree{0};
        std::size_t coeff_count_power{0};
        uint64_t modulus{0};

        /*
        * Storage used by the standalone constructor.
        *
        * These vectors remain empty when externally owned tables are used.
        */
        std::vector<NTTMultiplyOperand> owned_root_powers;
        std::vector<NTTMultiplyOperand> owned_inv_root_powers;

        /*
        * Active tables used by the transform kernels.
        *
        * They point either to the owned vectors above or to externally
        * supplied tables.
        */
        const NTTMultiplyOperand *root_powers{nullptr};
        const NTTMultiplyOperand *inv_root_powers{nullptr};

        NTTMultiplyOperand inv_degree_modulo{0, 0};
    };

    namespace {

        size_t ComputeCoeffCountPower(uint64_t degree) {
            size_t result = 0;
            while (degree > 1) {
                degree >>= 1;
                ++result;
            }
            return result;
        }

        /*
        * Internal forward-lazy dispatcher.
        */
        void NTTNegacyclicHarveyLazy(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers) {
            #ifdef AXHEL_HAS_SVE
            AXHEL_LOG("NTTNegacyclicHarveyLazy -> SVE" << ", degree=" << (size_t{1} << coeff_count_power));
            NTTNegacyclicHarveyLazySVE(operand, coeff_count_power, modulus, root_powers);
            #else
            AXHEL_LOG("NTTNegacyclicHarveyLazy -> native" << ", degree=" << (size_t{1} << coeff_count_power));
            NTTNegacyclicHarveyLazyNative(operand, coeff_count_power, modulus, root_powers);
            #endif
        }

        /*
        * Internal normalized forward NTT.
        */
        void NTTNegacyclicHarvey(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *root_powers) {
            NTTNegacyclicHarveyLazy(operand, coeff_count_power, modulus, root_powers);

            const size_t coeff_count = size_t{ 1 } << coeff_count_power;

            /*
            * Forward lazy output:
            *
            *   [0, 4q) -> [0, q)
            */
            EltwiseReduceMod(operand, operand, static_cast<uint64_t>(coeff_count), modulus, 4, 1);
        }

        /*
        * Internal inverse-lazy dispatcher.
        *
        */
        void InverseNTTNegacyclicHarveyLazy(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo) {
            #ifdef AXHEL_HAS_SVE
            AXHEL_LOG("InverseNTTNegacyclicHarveyLazy -> SVE" << ", degree=" << (size_t{1} << coeff_count_power));
            InverseNTTNegacyclicHarveyLazySVE(operand, coeff_count_power, modulus, inv_root_powers, inv_degree_modulo);
            #else
            AXHEL_LOG("InverseNTTNegacyclicHarveyLazy -> native" << ", degree=" << (size_t{1} << coeff_count_power));
            InverseNTTNegacyclicHarveyLazyNative(operand, coeff_count_power, modulus, inv_root_powers, inv_degree_modulo);
            #endif
        }

        /*
        * Internal normalized inverse NTT.
        */
        void InverseNTTNegacyclicHarvey(uint64_t *operand, size_t coeff_count_power, uint64_t modulus, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo) {
            InverseNTTNegacyclicHarveyLazy(operand, coeff_count_power, modulus, inv_root_powers, inv_degree_modulo);

            const size_t coeff_count = size_t{ 1 } << coeff_count_power;

            /*
            * Inverse lazy output:
            *
            *   [0, 2q) -> [0, q)
            */
            EltwiseReduceMod(operand, operand, static_cast<uint64_t>(coeff_count), modulus, 2, 1);
        }
    
    } // namespace anonymous 


    /// @brief Constructs an NTT instance from precomputed twiddle-factor tables.
    NTT::NTT(uint64_t degree, uint64_t modulus, const NTTMultiplyOperand *root_powers, const NTTMultiplyOperand *inv_root_powers, NTTMultiplyOperand inv_degree_modulo): impl_(std::make_shared<Impl>()) {
        impl_->degree = degree;
        impl_->coeff_count_power = ComputeCoeffCountPower(degree);
        impl_->modulus = modulus;
        impl_->root_powers = root_powers;
        impl_->inv_root_powers = inv_root_powers;
        impl_->inv_degree_modulo = inv_degree_modulo;

        AXHEL_LOG("NTT constructed with external tables" << ", degree=" << degree << ", modulus=" << modulus);
    }


    NTT::NTT(uint64_t degree, uint64_t modulus, uint64_t root_of_unity) : impl_(std::make_shared<Impl>()) {
        impl_->degree = degree;
        impl_->coeff_count_power = ComputeCoeffCountPower(degree);
        impl_->modulus = modulus;

        detail::NTTTableData tables = detail::GenerateNTTTables(static_cast<std::size_t>(degree), impl_->coeff_count_power, modulus, root_of_unity);

        impl_->owned_root_powers = std::move(tables.root_powers);
        impl_->owned_inv_root_powers = std::move(tables.inv_root_powers);
        impl_->inv_degree_modulo = tables.inv_degree_modulo;

        /*
        * Set active pointers only after the vectors have received their
        * final storage.
        */
        impl_->root_powers = impl_->owned_root_powers.data();
        impl_->inv_root_powers = impl_->owned_inv_root_powers.data();

        AXHEL_LOG("NTT constructed with internally generated tables" << ", degree=" << degree << ", modulus=" << modulus << ", root=" << root_of_unity);
    }

    
    NTT::~NTT() = default;


    NTT::NTT(const NTT &) noexcept = default;


    NTT &NTT::operator=(const NTT &) noexcept = default;


    NTT::NTT(NTT &&) noexcept = default;


    NTT &NTT::operator=(NTT &&) noexcept = default;

    /// @brief Computes the forward negacyclic Harvey NTT.
    void NTT::ComputeForward(uint64_t *result, const uint64_t *operand, uint64_t input_mod_factor, uint64_t output_mod_factor) const {
        if (result != operand) {
         std::copy_n(operand, static_cast<std::size_t>(impl_->degree), result);
        }

        if (output_mod_factor == 4) {
            NTTNegacyclicHarveyLazy(result, impl_->coeff_count_power, impl_->modulus, impl_->root_powers);
        }
        else {
            NTTNegacyclicHarvey(result, impl_->coeff_count_power, impl_->modulus, impl_->root_powers);
        }
    }


    /// @brief Computes the inverse negacyclic Harvey NTT.
    void NTT::ComputeInverse(uint64_t *result, const uint64_t *operand, uint64_t input_mod_factor, uint64_t output_mod_factor) const {
        if (result != operand) {
            std::copy_n(operand, static_cast<std::size_t>(impl_->degree), result);
        }
        if (output_mod_factor == 2) {
            InverseNTTNegacyclicHarveyLazy(result, impl_->coeff_count_power, impl_->modulus, impl_->inv_root_powers, impl_->inv_degree_modulo);
        }
        else {
            InverseNTTNegacyclicHarvey(result, impl_->coeff_count_power, impl_->modulus, impl_->inv_root_powers, impl_->inv_degree_modulo);
        }
    }


    /*
    * Accessors 
    */
    uint64_t NTT::Degree() const noexcept {
        return impl_ ? impl_->degree : 0;
    }

    uint64_t NTT::Modulus() const noexcept {
        return impl_ ? impl_->modulus : 0;
    }


} // namespace axhel
} // namespace unipi