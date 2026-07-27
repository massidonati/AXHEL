// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "ntt/ntt-sve.hpp"
#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/sve-arith.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "ntt/ntt-internal.hpp"
#include "axhel/util/compiler.hpp"
#include <stddef.h>
#include <stdint.h>

#ifdef AXHEL_HAS_SVE

#include <arm_sve.h>


namespace unipi {
namespace axhel {

namespace {

    inline void ForwardButterflyVector(
        svbool_t pg,
        svuint64_t x,
        svuint64_t y,
        svuint64_t root,
        svuint64_t root_quotient,
        svuint64_t modulus,
        svuint64_t twice_modulus,
        svuint64_t &result_x,
        svuint64_t &result_y) noexcept
    {

        const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg, x, twice_modulus);
        const svuint64_t transformed_y = MultiplyUIntModLazySVE(pg, y, root, root_quotient, modulus);
        result_x = svadd_u64_x(pg, guarded_x, transformed_y);
        result_y = svsub_u64_x(pg, svadd_u64_x(pg, guarded_x, twice_modulus), transformed_y);
    }


    inline void InverseButterflyVector(
        svbool_t pg,
        svuint64_t x,
        svuint64_t y,
        svuint64_t inv_root,
        svuint64_t inv_root_quotient,
        svuint64_t modulus,
        svuint64_t twice_modulus,
        svuint64_t &result_x,
        svuint64_t &result_y) noexcept
    {

        const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg, x, twice_modulus);
        result_x = svadd_u64_x(pg, guarded_x, y);
        const svuint64_t difference = svsub_u64_x(pg, svadd_u64_x(pg, guarded_x, twice_modulus), y);
        result_y = MultiplyUIntModLazySVE(pg, difference, inv_root, inv_root_quotient, modulus);
    }


    inline void ForwardStageScalar(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {

        for (size_t i = 0; i < m; ++i) {
            const NTTMultiplyOperand &root = root_powers[m + i];
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            AXHEL_UNROLL(4)
            for (size_t j = 0; j < gap; ++j) {
                detail::ForwardButterflyNative(x_ptr[j], y_ptr[j], root, modulus, twice_modulus);
            }
        }
    }

    inline void ForwardStageSVE(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {

        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);
        const size_t vectorized_count = gap - (gap % lanes);

        for (size_t i = 0; i < m; ++i) {
            const NTTMultiplyOperand &root = root_powers[m + i];
            const svuint64_t root_vec = svdup_n_u64(root.operand);
            const svuint64_t root_quotient_vec = svdup_n_u64(root.quotient);
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr =  x_ptr + gap;

            size_t j = 0;

            /*
            * Fully populated SVE vectors.
            */
            AXHEL_UNROLL(4)
            for (; j < vectorized_count; j += lanes) {
                const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                ForwardButterflyVector(pg_all, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_all, x_ptr + j, result_x);
                svst1_u64(pg_all, y_ptr + j, result_y);
            }

            /*
            * At most one partial vector.
            */
            if (j < gap) {
                const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(j), static_cast<uint64_t>(gap));
                const svuint64_t x = svld1_u64(pg_tail, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_tail, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                ForwardButterflyVector(pg_tail, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_tail, x_ptr + j, result_x);
                svst1_u64(pg_tail, y_ptr + j, result_y);
            }
        }
    }

    inline void InverseStageScalar(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        size_t &root_index) noexcept
    {

        for (size_t i = 0; i < m; ++i) {
            const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            AXHEL_UNROLL(4)
            for (size_t j = 0; j < gap; ++j) {
                detail::InverseButterflyNative(x_ptr[j], y_ptr[j], inv_root, modulus, twice_modulus);
            }
        }
    }

    inline void InverseStageSVE(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        size_t &root_index) noexcept
    {

        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);
        const size_t vectorized_count = gap - (gap % lanes);

        for (size_t i = 0; i < m; ++i) {
            const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];
            const svuint64_t inv_root_vec = svdup_n_u64(inv_root.operand);
            const svuint64_t inv_root_quotient_vec = svdup_n_u64(inv_root.quotient);
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            size_t j = 0;

            AXHEL_UNROLL(4)
            for (; j < vectorized_count; j += lanes) {
                const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                InverseButterflyVector(pg_all, x, y, inv_root_vec, inv_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_all, x_ptr + j, result_x);
                svst1_u64(pg_all, y_ptr + j, result_y);
            }

            if (j < gap) {
                const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(j), static_cast<uint64_t>(gap));

                const svuint64_t x = svld1_u64(pg_tail, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_tail, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                InverseButterflyVector(pg_tail, x, y, inv_root_vec, inv_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_tail, x_ptr + j, result_x);
                svst1_u64(pg_tail, y_ptr + j, result_y);
            }
        }
    }


    inline void InverseFinalStageSVE(
        uint64_t *__restrict operand,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand &scaled_root,
        const NTTMultiplyOperand &inv_degree_modulo) noexcept
    {
        
        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);
        const svuint64_t inv_degree_vec = svdup_n_u64(inv_degree_modulo.operand);
        const svuint64_t inv_degree_quotient_vec = svdup_n_u64(inv_degree_modulo.quotient);
        const svuint64_t scaled_root_vec = svdup_n_u64(scaled_root.operand);
        const svuint64_t scaled_root_quotient_vec = svdup_n_u64(scaled_root.quotient);
        uint64_t *x_ptr = operand;
        uint64_t *y_ptr = operand + gap;
        const size_t vectorized_count = gap - (gap % lanes);

        size_t j = 0;

        AXHEL_UNROLL(4)
        for (; j < vectorized_count; j += lanes) {
            const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_all, y_ptr + j);
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_all, x, twice_modulus_vec);
            const svuint64_t sum = svadd_u64_x(pg_all, guarded_x, y);
            const svuint64_t difference = svsub_u64_x(pg_all, svadd_u64_x(pg_all, guarded_x, twice_modulus_vec), y);
            const svuint64_t result_x = MultiplyUIntModLazySVE(pg_all, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);
            const svuint64_t result_y = MultiplyUIntModLazySVE(pg_all, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            svst1_u64(pg_all, x_ptr + j, result_x);
            svst1_u64(pg_all, y_ptr + j, result_y);
        }

        if (j < gap) {
            const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(j), static_cast<uint64_t>(gap));
            const svuint64_t x = svld1_u64(pg_tail, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_tail, y_ptr + j);
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_tail, x, twice_modulus_vec); 
            const svuint64_t sum = svadd_u64_x(pg_tail, guarded_x, y);
            const svuint64_t difference = svsub_u64_x( pg_tail, svadd_u64_x(pg_tail, guarded_x, twice_modulus_vec), y);
            const svuint64_t result_x = MultiplyUIntModLazySVE(pg_tail, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);
            const svuint64_t result_y = MultiplyUIntModLazySVE(pg_tail, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            svst1_u64(pg_tail, x_ptr + j, result_x);
            svst1_u64(pg_tail, y_ptr + j, result_y);
        }
    }

} // namespace


    void NTTNegacyclicHarveyLazySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *root_powers)
    {
       
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        const size_t lanes = static_cast<size_t>(svcntd());

        size_t m = 1;
        size_t gap = coeff_count >> 1;

        while (m < coeff_count) {
            /*
            * Use SVE only when a block contains at least one complete
            * vector of independent butterflies.
            */
            if (gap >= lanes) {
                ForwardStageSVE(operand, m, gap, modulus, twice_modulus, root_powers);
            }
            else {
                ForwardStageScalar(operand, m, gap, modulus, twice_modulus, root_powers);
            }

            m <<= 1;
            gap >>= 1;
        }
    }

    
    void InverseNTTNegacyclicHarveyLazySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo)
    {

        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        const size_t lanes = static_cast<size_t>(svcntd());

        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        while (m > 1) {
            if (gap >= lanes) {
                InverseStageSVE(operand, m, gap, modulus, twice_modulus, inv_root_powers, root_index);
            } 
            else {
                InverseStageScalar(operand, m, gap, modulus, twice_modulus, inv_root_powers, root_index);
            }

            m >>= 1;
            gap <<= 1;
        }

        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{
            scaled_root_operand,
            ComputeShoupQuotient(scaled_root_operand, modulus)
        };

        InverseFinalStageSVE(operand, gap, modulus, twice_modulus, scaled_root, inv_degree_modulo);
    }

} // namespace axhel
} // namespace unipi

#endif