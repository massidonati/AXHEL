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
        
        const svuint64_t sum = svadd_u64_x(pg, x, y);
        result_x = ReduceModFactor4To2SVE(pg, sum, twice_modulus);
        const svuint64_t difference = svsub_u64_x(pg, svadd_u64_x(pg, x, twice_modulus), y);
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


    inline void ForwardFinalStageSVE(
        uint64_t *__restrict operand,
        size_t m,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * Final forward stage:
        *
        *   gap = 1
        *
        * Memory layout:
        *
        *   x0, y0, x1, y1, ...
        *
        * Each butterfly uses a different root:
        *
        *   root_powers[m + i]
        */
        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);
        size_t i = 0;

        /*
        * Each SVE vector processes `lanes` independent butterflies.
        *
        * svld2_u64 works for every supported SVE vector length:
        *
        *   VL=128 -> 2 butterflies
        *   VL=256 -> 4 butterflies
        *   VL=512 -> 8 butterflies
        */
        AXHEL_UNROLL(4)
        for (; i + lanes <= m; i += lanes) {
            const svuint64x2_t xy = svld2_u64(pg_all, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            /*
            * NTTMultiplyOperand contains:
            *
            *   operand
            *   quotient
            *
            * Consecutive entries therefore form an interleaved sequence:
            *
            *   root0, quotient0, root1, quotient1, ...
            */
            const uint64_t *root_words = reinterpret_cast<const uint64_t *>(root_powers + m + i);
            const svuint64x2_t roots = svld2_u64(pg_all, root_words);
            const svuint64_t root_vec = svget2_u64(roots, 0);
            const svuint64_t root_quotient_vec = svget2_u64(roots, 1);

            svuint64_t result_x;
            svuint64_t result_y;

            ForwardButterflyVector(pg_all, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

            const svuint64x2_t results = svcreate2_u64(result_x, result_y);

            svst2_u64(pg_all, operand + (i << 1), results);
        }

        /*
        * Generic tail.
        *
        * For the usual power-of-two NTT sizes 
        * this branch should not be entered because m is divisible
        * by svcntd().
        */
        if (i < m) {
            const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(i), static_cast<uint64_t>(m));
            const svuint64x2_t xy = svld2_u64(pg_tail, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);
            const uint64_t *root_words = reinterpret_cast<const uint64_t *>(root_powers + m + i);
            const svuint64x2_t roots = svld2_u64(pg_tail, root_words);
            const svuint64_t root_vec = svget2_u64(roots, 0);
            const svuint64_t root_quotient_vec = svget2_u64(roots, 1);

            svuint64_t result_x;
            svuint64_t result_y;

            ForwardButterflyVector(pg_tail, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

            const svuint64x2_t results =svcreate2_u64(result_x, result_y);

            svst2_u64(pg_tail, operand + (i << 1), results);
        }
    }


    inline void ForwardFinalTwoStagesSVE128(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * This kernel fuses the final two forward stages:
        *
        *   first stage:  gap = 2, m = coeff_count / 4
        *   final stage:  gap = 1, m = coeff_count / 2
        *
        * It is specialized for a 128-bit SVE vector length:
        *
        *   svcntd() == 2
        *
        * Each iteration processes four consecutive coefficients:
        *
        *   a, b, c, d
        *
        * The gap = 2 stage applies butterflies to:
        *
        *   (a, c)
        *   (b, d)
        *
        * The gap = 1 stage then applies butterflies to:
        *
        *   (u0, u1)
        *   (v0, v1)
        *
        * where:
        *
        *   [u0, u1] and [v0, v1]
        *
        * are the results of the first fused stage.
        */
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        /*
        * Twiddle-table offsets for the two fused stages.
        */
        const size_t first_stage_m = coeff_count >> 2;
        const size_t final_stage_m = coeff_count >> 1;

        AXHEL_UNROLL(4)
        for (size_t i = 0; i < first_stage_m; ++i) {
            uint64_t *block = operand + (i << 2);

            /*
            * Four contiguous coefficients:
            *
            *   block[0] = a
            *   block[1] = b
            *   block[2] = c
            *   block[3] = d
            *
            * With VL=128:
            *
            *   x = [a, b]
            *   y = [c, d]
            *
            * These are exactly the two independent butterflies
            * of the gap = 2 stage.
            */
            const svuint64_t x = svld1_u64(pg_all, block);
            const svuint64_t y = svld1_u64(pg_all, block + 2);

            /*
            * Both lanes use the same twiddle in the gap = 2 stage.
            */
            const NTTMultiplyOperand &first_root = root_powers[first_stage_m + i];
            const svuint64_t first_root_vec = svdup_n_u64(first_root.operand);
            const svuint64_t first_root_quotient_vec = svdup_n_u64(first_root.quotient);

            svuint64_t upper;
            svuint64_t lower;

            ForwardButterflyVector(pg_all, x, y, first_root_vec, first_root_quotient_vec, modulus_vec, twice_modulus_vec, upper, lower);

            /*
            * After the gap = 2 stage:
            *
            *   upper = [u0, u1]
            *   lower = [v0, v1]
            *
            * The final gap = 1 stage needs the pairs:
            *
            *   (u0, u1)
            *   (v0, v1)
            *
            * Arrange them by transposing the two vectors:
            *
            *   final_x = [u0, v0]
            *   final_y = [u1, v1]
            *
            * Therefore:
            *
            *   lane 0 executes (u0, u1)
            *   lane 1 executes (v0, v1)
            */
            const svuint64_t final_x = svtrn1_u64(upper, lower);
            const svuint64_t final_y = svtrn2_u64(upper, lower);

            /*
            * The two final butterflies use two different twiddles:
            *
            *   root_powers[final_stage_m + 2*i]
            *   root_powers[final_stage_m + 2*i + 1]
            *
            * NTTMultiplyOperand is laid out as:
            *
            *   operand, quotient
            *
            * Thus svld2 deinterleaves the two structures into:
            *
            *   roots     = [root0, root1]
            *   quotients = [quotient0, quotient1]
            */
            const uint64_t *final_root_words = reinterpret_cast<const uint64_t *>(root_powers + final_stage_m + (i << 1));
            const svuint64x2_t final_roots = svld2_u64(pg_all, final_root_words);
            const svuint64_t final_root_vec = svget2_u64(final_roots, 0);
            const svuint64_t final_root_quotient_vec = svget2_u64( final_roots, 1);

            svuint64_t result_x;
            svuint64_t result_y;

            ForwardButterflyVector(pg_all, final_x, final_y, final_root_vec, final_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

            /*
            * The second butterfly produces:
            *
            *   result_x = [r0, r2]
            *   result_y = [r1, r3]
            *
            * Restore the contiguous order:
            *
            *   r0, r1, r2, r3
            */
            const svuint64_t output01 = svzip1_u64(result_x, result_y);
            const svuint64_t output23 = svzip2_u64(result_x, result_y);
           
            svst1_u64(pg_all, block, output01);
            svst1_u64(pg_all, block + 2, output23);
        }
    }


    /*
     * Generic forward NTT implementation for a compile-time
     * number of uint64_t SVE lanes.
     *
     * Lanes == 0 means that the vector length must be obtained
     * at runtime through svcntd().
     *
     * The generic path processes all stages down to gap = 1,
     * then uses the vector-length-independent final-stage kernel.
     */
    template <size_t Lanes>
    inline void ForwardTransformSVE(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        
        const size_t lanes = Lanes != 0 ? Lanes : static_cast<size_t>(svcntd());

        size_t m = 1;
        size_t gap = coeff_count >> 1;

        /*
         * Process every stage except the final gap = 1 stage.
         */
        while (gap > 1) {
            if (gap >= lanes) {
                ForwardStageSVE(operand, m, gap, modulus, twice_modulus, root_powers);
            }
            else {
                /*
                 * This is required only when the number of
                 * coefficients in a butterfly block is smaller
                 * than the current SVE vector length.
                 */
                ForwardStageScalar(operand, m, gap, modulus, twice_modulus, root_powers);
            }

            m <<= 1;
            gap >>= 1;
        }

        /*
         * Vector-length-independent final stage introduced
         */
        ForwardFinalStageSVE(operand, m, modulus, twice_modulus, root_powers);
    }

    
    /*
     * Specialized forward NTT implementation for:
     *
     *   SVE VL = 128 bits
     *   svcntd() = 2
     *
     * The final gap = 2 and gap = 1 stages are fused by
     * ForwardFinalTwoStagesSVE128.
     */
    template <>
    inline void ForwardTransformSVE<2>(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        
        size_t m = 1;
        size_t gap = coeff_count >> 1;

        /*
         * Process all stages preceding:
         *
         *   gap = 2
         *   gap = 1
         */
        while (gap > 2) {
            ForwardStageSVE(operand, m, gap, modulus, twice_modulus, root_powers);

            m <<= 1;
            gap >>= 1;
        }

        /*
         * A transform with at least four coefficients reaches
         * the fused final-two-stage kernel.
         */
        if (coeff_count >= 4) {
            ForwardFinalTwoStagesSVE128(operand, coeff_count, modulus, twice_modulus, root_powers);
        }
        else {
            /*
             * Defensive fallback for very small transforms.
             * Normal SEAL parameter sets never enter this branch.
             */
            ForwardFinalStageSVE(operand, m, modulus, twice_modulus, root_powers);
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

        /*
        * Runtime selection of a compile-time-specialized kernel.
        */
        switch (lanes) {
        case 2:
            /*
            * SVE VL=128:
            */
            ForwardTransformSVE<2>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;

        default:
            /*
            * Vector-length-agnostic fallback for all other SVE
            * vector lengths.
            */
            ForwardTransformSVE<0>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;
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