// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include "ntt/ntt-sve.hpp"
#include "axhel/eltwise/eltwise-reduce-mod.hpp"
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

    /* ********************************* */
    /*  FORWARD NTT                      */
    /* ********************************* */

    /// @brief Performs a lane-wise forward Harvey butterfly using SVE.
    /// @pre Each active lane of x and y must be in the range [0, 4q).
    /// @post Each active lane of result_x and result_y is in the range [0, 4q).
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
        // [0, 4q) -> [0, 2q).
        const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg, x, twice_modulus);

         // Lazy Shoup product in [0, 2q).
        const svuint64_t transformed_y = MultiplyUIntModLazySVE(pg, y, root, root_quotient, modulus);

        // Results remain in [0, 4q).
        result_x = svadd_u64_x(pg, guarded_x, transformed_y);
        result_y = svsub_u64_x(pg, svadd_u64_x(pg, guarded_x, twice_modulus), transformed_y);
    }


    /// @brief Computes one forward NTT stage using scalar Harvey butterflies.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    inline void ForwardStageScalar(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {

        for (size_t i = 0; i < m; ++i) {
            // Forward roots are stored in bit-reversed order.
            const NTTMultiplyOperand &root = root_powers[m + i];

            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            AXHEL_UNROLL(4)
            // Apply the forward Harvey butterflies for the current block.
            for (size_t j = 0; j < gap; ++j) {
                detail::ForwardButterflyNative(x_ptr[j], y_ptr[j], root, modulus, twice_modulus);
            }
        }
    }


    /// @brief Computes one forward NTT stage using SVE Harvey butterflies.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
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
             // Forward roots are stored in bit-reversed order.
            const NTTMultiplyOperand &root = root_powers[m + i];

            const svuint64_t root_vec = svdup_n_u64(root.operand);
            const svuint64_t root_quotient_vec = svdup_n_u64(root.quotient);
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr =  x_ptr + gap;

            size_t j = 0;

            AXHEL_UNROLL(4)
            // Process full-width SVE vectors.
            for (; j < vectorized_count; j += lanes) {
                const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                ForwardButterflyVector(pg_all, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_all, x_ptr + j, result_x);
                svst1_u64(pg_all, y_ptr + j, result_y);
            }

             // Handle the remaining elements with predication.
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


    /// @brief Computes one forward NTT stage specialized for SVE VL=128.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    inline void ForwardStageSVE128(
        uint64_t *__restrict operand,
        size_t m,
        size_t gap,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        
        constexpr size_t lanes = 2;
        constexpr size_t ilp_step = 2 * lanes;
        
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        for (size_t i = 0; i < m; ++i) {
             // Forward roots are stored in bit-reversed order.
            const NTTMultiplyOperand &root = root_powers[m + i];

            const svuint64_t root_vec = svdup_n_u64(root.operand);
            const svuint64_t root_quotient_vec =  svdup_n_u64(root.quotient);
            
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            size_t j = 0;

            AXHEL_UNROLL(2)
            // Process two full-width SVE vectors per iteration.
            for (; j + ilp_step <= gap; j += ilp_step) {
                const svuint64_t x0 = svld1_u64(pg_all, x_ptr + j);
                const svuint64_t y0 = svld1_u64(pg_all, y_ptr + j);
                const svuint64_t x1 = svld1_u64(pg_all, x_ptr + j + lanes);
                const svuint64_t y1 = svld1_u64(pg_all, y_ptr + j + lanes);

                svuint64_t result_x0;
                svuint64_t result_y0;
                svuint64_t result_x1;
                svuint64_t result_y1;

                ForwardButterflyVector(pg_all, x0, y0, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x0, result_y0);
                ForwardButterflyVector( pg_all, x1, y1, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x1, result_y1);

                svst1_u64(pg_all, x_ptr + j, result_x0);
                svst1_u64(pg_all, y_ptr + j, result_y0);
                svst1_u64(pg_all, x_ptr + j + lanes, result_x1);
                svst1_u64(pg_all, y_ptr + j + lanes, result_y1);
            }
 
            // Process one remaining full-width SVE vector.
            if (j + lanes <= gap) {
                const svuint64_t x = svld1_u64( pg_all, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                ForwardButterflyVector(pg_all, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_all, x_ptr + j, result_x);
                svst1_u64(pg_all, y_ptr + j, result_y);

                j += lanes;
            }

            // Handle any remaining elements with the scalar butterfly.
            for (; j < gap; ++j) {
                detail::ForwardButterflyNative(x_ptr[j], y_ptr[j], root, modulus, twice_modulus);
            }
        }
    }


    /// @brief Computes the final forward NTT stage using SVE.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    inline void ForwardFinalStageSVE(
        uint64_t *__restrict operand,
        size_t m,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * Final forward stage with gap = 1.
        *
        * Consecutive coefficients form butterfly pairs:
        *
        *   x0, y0, x1, y1, ...
        */
        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);
        size_t i = 0;

        AXHEL_UNROLL(4)
        // Process lanes independent butterflies per iteration.
        //  VL=128 -> 2 butterflies; VL=256 -> 4 butterflies; VL=512 -> 8 butterflies
        for (; i + lanes <= m; i += lanes) {
            const svuint64x2_t xy = svld2_u64(pg_all, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            /*
            * Load consecutive roots and Shoup quotients:
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

        // Handle the remaining butterflies with predication.
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


    /// @brief Computes the normalized final forward NTT stage using SVE.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, modulus).
    inline void ForwardFinalStageNormalizedSVE(
        uint64_t *__restrict operand,
        size_t m,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * Final forward stage with gap = 1.
        *
        * Consecutive coefficients form butterfly pairs:
        *
        *   x0, y0, x1, y1, ...
        */
        const size_t lanes = static_cast<size_t>(svcntd());
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        size_t i = 0;

        AXHEL_UNROLL(4)
        // Process lanes independent butterflies per iteration.
        for (; i + lanes <= m; i += lanes) {
            const svuint64x2_t xy = svld2_u64(pg_all, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            /*
            * Load consecutive roots and Shoup quotients:
            *
            *   root0, quotient0, root1, quotient1, ...
            */
            const uint64_t *root_words = reinterpret_cast<const uint64_t *>(root_powers + m + i);
            const svuint64x2_t roots = svld2_u64(pg_all, root_words);
            const svuint64_t root_vec = svget2_u64(roots, 0);
            const svuint64_t root_quotient_vec = svget2_u64(roots, 1);

            svuint64_t result_x_lazy;
            svuint64_t result_y_lazy;

            ForwardButterflyVector(pg_all, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x_lazy, result_y_lazy);

            // Normalize the lazy results to [0, q).
            const svuint64_t result_x_2q = ReduceModFactor4To2SVE(pg_all, result_x_lazy, twice_modulus_vec);
            const svuint64_t result_y_2q = ReduceModFactor4To2SVE(pg_all, result_y_lazy, twice_modulus_vec);
            const svuint64_t result_x = ReduceModFactor2To1SVE(pg_all, result_x_2q, modulus_vec);
            const svuint64_t result_y = ReduceModFactor2To1SVE(pg_all, result_y_2q, modulus_vec);

            const svuint64x2_t results = svcreate2_u64(result_x, result_y);

            svst2_u64(pg_all, operand + (i << 1), results);
        }

        // Handle the remaining butterflies with predication.
        if (i < m) {
            const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(i), static_cast<uint64_t>(m));
            const svuint64x2_t xy = svld2_u64(pg_tail, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            const uint64_t *root_words = reinterpret_cast<const uint64_t *>(root_powers + m + i);
            const svuint64x2_t roots = svld2_u64(pg_tail, root_words);
            const svuint64_t root_vec = svget2_u64(roots, 0);
            const svuint64_t root_quotient_vec = svget2_u64(roots, 1);

            svuint64_t result_x_lazy;
            svuint64_t result_y_lazy;

            ForwardButterflyVector(pg_tail, x, y, root_vec, root_quotient_vec, modulus_vec, twice_modulus_vec, result_x_lazy, result_y_lazy);

            // Normalize the lazy results to [0, q).
            const svuint64_t result_x_2q = ReduceModFactor4To2SVE(pg_tail, result_x_lazy, twice_modulus_vec);
            const svuint64_t result_y_2q = ReduceModFactor4To2SVE(pg_tail, result_y_lazy, twice_modulus_vec);
            const svuint64_t result_x = ReduceModFactor2To1SVE(pg_tail, result_x_2q, modulus_vec);
            const svuint64_t result_y = ReduceModFactor2To1SVE(pg_tail, result_y_2q, modulus_vec);

            const svuint64x2_t results = svcreate2_u64(result_x, result_y);

            svst2_u64(pg_tail, operand + (i << 1), results);
        }
    }


    /// @brief Computes the final two forward NTT stages specialized for SVE VL=128.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    inline void ForwardFinalTwoStagesSVE128(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * Fuse the final two forward stages:
        *
        *   gap = 2, m = coeff_count / 4
        *   gap = 1, m = coeff_count / 2
        *
        * VL=128 provides two 64-bit lanes.
        */
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        // Twiddle-table offsets for the two fused stages.
        const size_t first_stage_m = coeff_count >> 2;
        const size_t final_stage_m = coeff_count >> 1;

        AXHEL_UNROLL(4)
        // Process two gap = 2 butterflies per iteration.
        for (size_t i = 0; i < first_stage_m; ++i) {
            uint64_t *block = operand + (i << 2);

            /*
            * Load four coefficients as:
            *
            *   x = [a, b]
            *   y = [c, d]
            *
            * corresponding to the two gap = 2 butterflies.
            */
            const svuint64_t x = svld1_u64(pg_all, block);
            const svuint64_t y = svld1_u64(pg_all, block + 2);

            // Both lanes use the same root in the gap = 2 stage.
            const NTTMultiplyOperand &first_root = root_powers[first_stage_m + i];
            const svuint64_t first_root_vec = svdup_n_u64(first_root.operand);
            const svuint64_t first_root_quotient_vec = svdup_n_u64(first_root.quotient);

            svuint64_t upper;
            svuint64_t lower;

            ForwardButterflyVector(pg_all, x, y, first_root_vec, first_root_quotient_vec, modulus_vec, twice_modulus_vec, upper, lower);

            /*
            * Rearrange the first-stage results for the final gap = 1 stage:
            *
            *   upper = [u0, u1]
            *   lower = [v0, v1]
            *
            *   final_x = [u0, v0]
            *   final_y = [u1, v1]
            */
            const svuint64_t final_x = svtrn1_u64(upper, lower);
            const svuint64_t final_y = svtrn2_u64(upper, lower);

            // Load the two roots required by the final-stage butterflies.
            const uint64_t *final_root_words = reinterpret_cast<const uint64_t *>(root_powers + final_stage_m + (i << 1));
            const svuint64x2_t final_roots = svld2_u64(pg_all, final_root_words);
            const svuint64_t final_root_vec = svget2_u64(final_roots, 0);
            const svuint64_t final_root_quotient_vec = svget2_u64( final_roots, 1);

            svuint64_t result_x;
            svuint64_t result_y;

            ForwardButterflyVector(pg_all, final_x, final_y, final_root_vec, final_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);
           
            /*
            * Restore contiguous coefficient order:
            *
            *   result_x = [r0, r2]
            *   result_y = [r1, r3]
            *
            *   output = [r0, r1, r2, r3]
            */
            const svuint64_t output01 = svzip1_u64(result_x, result_y);
            const svuint64_t output23 = svzip2_u64(result_x, result_y);
           
            svst1_u64(pg_all, block, output01);
            svst1_u64(pg_all, block + 2, output23);
        }
    }


    /// @brief Computes and normalizes the final two forward NTT stages specialized for SVE VL=128.
    /// @pre Each input coefficient must be in the range [0, 4 * modulus).
    /// @post Each output coefficient is in the range [0, modulus).
    inline void ForwardFinalTwoStagesNormalizedSVE128(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        /*
        * Fuse the final two forward stages:
        *
        *   gap = 2, m = coeff_count / 4
        *   gap = 1, m = coeff_count / 2
        *
        * VL=128 provides two 64-bit lanes.
        */
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        const size_t first_stage_m = coeff_count >> 2;
        const size_t final_stage_m = coeff_count >> 1;

        AXHEL_UNROLL(4)
        // Process two gap = 2 butterflies per iteration.
        for (size_t i = 0; i < first_stage_m; ++i) {
            uint64_t *block = operand + (i << 2);

            /*
            * Load four coefficients as:
            *
            *   x = [a, b]
            *   y = [c, d]
            *
            * corresponding to the two gap = 2 butterflies.
            */
            const svuint64_t x = svld1_u64(pg_all, block);
            const svuint64_t y = svld1_u64(pg_all, block + 2);

            // Both lanes use the same root in the gap = 2 stage.
            const NTTMultiplyOperand &first_root = root_powers[first_stage_m + i];
            const svuint64_t first_root_vec = svdup_n_u64(first_root.operand);
            const svuint64_t first_root_quotient_vec = svdup_n_u64(first_root.quotient);

            svuint64_t upper;
            svuint64_t lower;

            ForwardButterflyVector(pg_all, x, y, first_root_vec, first_root_quotient_vec, modulus_vec, twice_modulus_vec, upper, lower);

            /*
            * Rearrange the first-stage results for the final gap = 1 stage:
            *
            *   upper = [u0, u1]
            *   lower = [v0, v1]
            *
            *   final_x = [u0, v0]
            *   final_y = [u1, v1]
            */
            const svuint64_t final_x = svtrn1_u64(upper, lower);
            const svuint64_t final_y = svtrn2_u64(upper, lower);

            // Load the two roots required by the final-stage butterflies.
            const uint64_t *final_root_words = reinterpret_cast<const uint64_t *>(root_powers + final_stage_m + (i << 1));
            const svuint64x2_t final_roots = svld2_u64(pg_all, final_root_words);
            const svuint64_t final_root_vec = svget2_u64(final_roots, 0);
            const svuint64_t final_root_quotient_vec = svget2_u64(final_roots, 1);

            svuint64_t result_x_lazy;
            svuint64_t result_y_lazy;

            ForwardButterflyVector(pg_all, final_x, final_y, final_root_vec, final_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x_lazy, result_y_lazy);

            // Normalize the lazy results to [0, q).
            const svuint64_t result_x_2q = ReduceModFactor4To2SVE(pg_all, result_x_lazy, twice_modulus_vec);
            const svuint64_t result_y_2q = ReduceModFactor4To2SVE(pg_all, result_y_lazy, twice_modulus_vec);
            const svuint64_t result_x = ReduceModFactor2To1SVE(pg_all, result_x_2q, modulus_vec);
            const svuint64_t result_y = ReduceModFactor2To1SVE(pg_all, result_y_2q, modulus_vec);

            /*
            * Restore contiguous coefficient order:
            *
            *   result_x = [r0, r2]
            *   result_y = [r1, r3]
            *
            *   output = [r0, r1, r2, r3]
            */
            const svuint64_t output01 = svzip1_u64(result_x, result_y);
            const svuint64_t output23 = svzip2_u64(result_x, result_y);

            svst1_u64(pg_all, block, output01);
            svst1_u64(pg_all, block + 2, output23);
        }
    }


    /// @brief Computes the forward NTT using SVE with a compile-time or runtime vector length.
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

        // Process all stages except the final gap = 1 stage.
        while (gap > 1) {
            if (gap >= lanes) {
                ForwardStageSVE(operand, m, gap, modulus, twice_modulus, root_powers);
            }
            else {
                // Use the scalar stage when the gap is smaller than the SVE vector length.
                ForwardStageScalar(operand, m, gap, modulus, twice_modulus, root_powers);
            }

            m <<= 1;
            gap >>= 1;
        }

        // Process the final gap = 1 stage.
        ForwardFinalStageSVE(operand, m, modulus, twice_modulus, root_powers);
    }

    
    /// @brief Computes the forward NTT specialized for SVE VL=128.
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

       // Process all stages preceding the final gap = 2 and gap = 1 stages.
        while (gap > 2) {
            ForwardStageSVE128(operand, m, gap, modulus, twice_modulus, root_powers);

            m <<= 1;
            gap >>= 1;
        }

        // Process the final two stages with the fused VL=128 kernel.
        if (coeff_count >= 4) {
            ForwardFinalTwoStagesSVE128(operand, coeff_count, modulus, twice_modulus, root_powers);
        }
        else {
            // Handle very small transforms with the generic final-stage kernel.
            ForwardFinalStageSVE(operand, m, modulus, twice_modulus, root_powers);
        }
    }


    /// @brief Computes the normalized forward NTT using SVE with a compile-time or runtime vector length.
    /// @note Lanes = 0 selects the runtime SVE vector length.
    template <size_t Lanes>
    inline void ForwardTransformNormalizedSVE(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        const size_t lanes = Lanes != 0 ? Lanes : static_cast<size_t>(svcntd());

        size_t m = 1;
        size_t gap = coeff_count >> 1;

        // Process all stages except the final gap = 1 stage.
        while (gap > 1) {
            if (gap >= lanes) {
                ForwardStageSVE(operand, m, gap, modulus, twice_modulus, root_powers);
            }
            else {
                // Use the scalar stage when the gap is smaller than the SVE vector length.
                ForwardStageScalar(operand, m, gap, modulus, twice_modulus, root_powers);
            }

            m <<= 1;
            gap >>= 1;
        }

        // Process and normalize the final gap = 1 stage.
        ForwardFinalStageNormalizedSVE(operand, m, modulus, twice_modulus, root_powers);
    }


    /// @brief Computes the normalized forward NTT specialized for SVE VL=128.
    /// @note Uses two 64-bit SVE lanes and fuses the final gap = 2 and gap = 1 stages.
    template <>
    inline void ForwardTransformNormalizedSVE<2>(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict root_powers) noexcept
    {
        size_t m = 1;
        size_t gap = coeff_count >> 1;

        // Process all stages preceding the final gap = 2 and gap = 1 stages.
        while (gap > 2) {
            ForwardStageSVE128(operand, m, gap, modulus, twice_modulus, root_powers);

            m <<= 1;
            gap >>= 1;
        }

        // Process and normalize the final two stages with the fused VL=128 kernel.
        if (coeff_count >= 4) {
            ForwardFinalTwoStagesNormalizedSVE128(operand, coeff_count, modulus, twice_modulus, root_powers);
        }
        else {
            // Handle very small transforms with the generic normalized final-stage kernel.
            ForwardFinalStageNormalizedSVE(operand, m, modulus, twice_modulus, root_powers);
        }
    }


    /* ********************************* */
    /*  INVERSE NTT                      */
    /* ********************************* */

    /// @brief Performs a lane-wise inverse Harvey butterfly using SVE.
    /// @pre Each active lane of x and y must be in the range [0, 2q).
    /// @post Each active lane of result_x and result_y is in the range [0, 2q).
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
        
        // [0, 2q) + [0, 2q) -> [0, 4q) -> [0, 2q).
        const svuint64_t sum = svadd_u64_x(pg, x, y);
        result_x = ReduceModFactor4To2SVE(pg, sum, twice_modulus);
        const svuint64_t difference = svsub_u64_x(pg, svadd_u64_x(pg, x, twice_modulus), y);

         // Lazy Shoup product in [0, 2q).
        result_y = MultiplyUIntModLazySVE(pg, difference, inv_root, inv_root_quotient, modulus);
    }


    /// @brief Computes one inverse NTT stage using scalar Harvey butterflies.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).
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
             // Inverse roots are consumed sequentially.
            const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];

            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            AXHEL_UNROLL(4)
             // Apply the inverse Harvey butterflies for the current block.
            for (size_t j = 0; j < gap; ++j) {
                detail::InverseButterflyNative(x_ptr[j], y_ptr[j], inv_root, modulus, twice_modulus);
            }
        }
    }


    /// @brief Computes one inverse NTT stage using SVE Harvey butterflies.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).    
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
            // Inverse roots are consumed sequentially.
            const NTTMultiplyOperand &inv_root = inv_root_powers[root_index++];
            const svuint64_t inv_root_vec = svdup_n_u64(inv_root.operand);
            const svuint64_t inv_root_quotient_vec = svdup_n_u64(inv_root.quotient);
            const size_t block_start = i * (gap << 1);
            uint64_t *x_ptr = operand + block_start;
            uint64_t *y_ptr = x_ptr + gap;

            size_t j = 0;

            AXHEL_UNROLL(4)
            // Process full-width SVE vectors.
            for (; j < vectorized_count; j += lanes) {
                const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
                const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

                svuint64_t result_x;
                svuint64_t result_y;

                InverseButterflyVector(pg_all, x, y, inv_root_vec, inv_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

                svst1_u64(pg_all, x_ptr + j, result_x);
                svst1_u64(pg_all, y_ptr + j, result_y);
            }

             // Handle the remaining elements with predication.
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


    /// @brief Computes the initial inverse NTT stage specialized for SVE VL=128.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).
    inline void InverseInitialStageSVE128(
        uint64_t *__restrict operand,
        size_t m,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        size_t &root_index) noexcept
    {
        /*
        * Initial inverse stage with gap = 1.
        *
        * Consecutive coefficients form butterfly pairs:
        *
        *   x0, y0, x1, y1, ...
        */
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        AXHEL_UNROLL(4)
        // Process two butterflies per iteration.
        for (size_t i = 0; i < m; i += 2) {

            /*
            * Load two butterfly pairs as:
            *
            *   x = [x0, x1]
            *   y = [y0, y1]
            */
            const svuint64x2_t xy = svld2_u64(pg_all, operand + (i << 1));
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            // Load two consecutive inverse roots and Shoup quotients.
            const uint64_t *root_words = reinterpret_cast<const uint64_t *>(inv_root_powers + root_index + i);
            const svuint64x2_t roots = svld2_u64(pg_all, root_words);
            const svuint64_t inv_root_vec = svget2_u64(roots, 0);
            const svuint64_t inv_root_quotient_vec = svget2_u64(roots, 1);

            svuint64_t result_x;
            svuint64_t result_y;

            InverseButterflyVector(pg_all, x, y, inv_root_vec, inv_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

            // Restore the interleaved coefficient layout.
            const svuint64x2_t results = svcreate2_u64(result_x, result_y);

            svst2_u64(pg_all, operand + (i << 1), results);
        }

        root_index += m;
    }


    /// @brief Computes the first two inverse NTT stages specialized for SVE VL=128.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).
    inline void InverseInitialTwoStagesSVE128(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        size_t &root_index) noexcept
    {
        /*
        * Fuse the first two inverse stages:
        *
        *   gap = 1, m = coeff_count / 2
        *   gap = 2, m = coeff_count / 4
        *
        * VL=128 provides two 64-bit lanes.
        */
        const svbool_t pg_all = svptrue_b64();
        const svuint64_t modulus_vec = svdup_n_u64(modulus);
        const svuint64_t twice_modulus_vec = svdup_n_u64(twice_modulus);

        const size_t first_stage_m = coeff_count >> 1;
        const size_t second_stage_m = coeff_count >> 2;

        // Root-table offsets for the two fused stages.
        const size_t first_stage_root_index = root_index;
        const size_t second_stage_root_index = root_index + first_stage_m;

        AXHEL_UNROLL(4)
        // Process one four-coefficient block per iteration.
        for (size_t i = 0; i < second_stage_m; ++i) {
            
            uint64_t *block = operand + (i << 2);

            /*
            * Load four coefficients as:
            *
            *   x = [a, c]
            *   y = [b, d]
            *
            * corresponding to the two gap = 1 butterflies.
            */
            const svuint64x2_t xy = svld2_u64(pg_all, block);
            const svuint64_t x = svget2_u64(xy, 0);
            const svuint64_t y = svget2_u64(xy, 1);

            // Load the two inverse roots required by the gap = 1 stage.
            const uint64_t *first_root_words = reinterpret_cast<const uint64_t *>(inv_root_powers + first_stage_root_index + (i << 1));
            const svuint64x2_t first_roots = svld2_u64(pg_all, first_root_words);
            const svuint64_t first_root_vec = svget2_u64(first_roots, 0);
            const svuint64_t first_root_quotient_vec = svget2_u64(first_roots, 1);

            svuint64_t first_x;
            svuint64_t first_y;

            InverseButterflyVector(pg_all, x, y, first_root_vec, first_root_quotient_vec, modulus_vec, twice_modulus_vec, first_x, first_y);

            /*
            * Rearrange the first-stage results for the gap = 2 stage:
            *
            *   first_x = [u0, u1]
            *   first_y = [v0, v1]
            *
            *   second_x = [u0, v0]
            *   second_y = [u1, v1]
            */
            const svuint64_t second_x = svtrn1_u64(first_x, first_y);
            const svuint64_t second_y = svtrn2_u64(first_x, first_y);

            // Both lanes use the same root in the gap = 2 stage.
            const NTTMultiplyOperand &second_root = inv_root_powers[second_stage_root_index + i];
            const svuint64_t second_root_vec = svdup_n_u64(second_root.operand);
            const svuint64_t second_root_quotient_vec = svdup_n_u64(second_root.quotient);

            svuint64_t result_x;
            svuint64_t result_y;

            InverseButterflyVector(pg_all, second_x, second_y, second_root_vec, second_root_quotient_vec, modulus_vec, twice_modulus_vec, result_x, result_y);

            // Store the gap = 2 stage results in contiguous order.
            svst1_u64(pg_all, block, result_x);
            svst1_u64(pg_all, block + 2, result_y);
        }
 
        // Advance past the roots consumed by both fused stages.
        root_index += first_stage_m + second_stage_m;
    }


    /// @brief Computes the final inverse NTT stage using SVE.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, 2 * modulus).
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
         // Process full-width SVE vectors.
        for (; j < vectorized_count; j += lanes) {
            const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

             // Guard x to the [0, 2q) range.
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_all, x, twice_modulus_vec);
            const svuint64_t sum = svadd_u64_x(pg_all, guarded_x, y);
            const svuint64_t difference = svsub_u64_x(pg_all, svadd_u64_x(pg_all, guarded_x, twice_modulus_vec), y);

            // Apply the inverse-degree scaling to the upper branch.
            const svuint64_t result_x = MultiplyUIntModLazySVE(pg_all, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);

             // Apply the scaled inverse root to the lower branch.
            const svuint64_t result_y = MultiplyUIntModLazySVE(pg_all, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            svst1_u64(pg_all, x_ptr + j, result_x);
            svst1_u64(pg_all, y_ptr + j, result_y);
        }

        // Handle the remaining elements with predication.
        if (j < gap) {
            const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(j), static_cast<uint64_t>(gap));
            const svuint64_t x = svld1_u64(pg_tail, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_tail, y_ptr + j);

             // Guard x to the [0, 2q) range.
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_tail, x, twice_modulus_vec); 
            const svuint64_t sum = svadd_u64_x(pg_tail, guarded_x, y);
            const svuint64_t difference = svsub_u64_x( pg_tail, svadd_u64_x(pg_tail, guarded_x, twice_modulus_vec), y);

            // Apply the inverse-degree scaling to the upper branch.
            const svuint64_t result_x = MultiplyUIntModLazySVE(pg_tail, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);

             // Apply the scaled inverse root to the lower branch.
            const svuint64_t result_y = MultiplyUIntModLazySVE(pg_tail, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            svst1_u64(pg_tail, x_ptr + j, result_x);
            svst1_u64(pg_tail, y_ptr + j, result_y);
        }
    }


    /// @brief Computes the normalized final inverse NTT stage using SVE.
    /// @pre Each input coefficient must be in the range [0, 2 * modulus).
    /// @post Each output coefficient is in the range [0, modulus).
    inline void InverseFinalStageNormalizedSVE(
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
         // Process full-width SVE vectors.
        for (; j < vectorized_count; j += lanes) {
            const svuint64_t x = svld1_u64(pg_all, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_all, y_ptr + j);

            // Guard x to the [0, 2q) range.
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_all, x, twice_modulus_vec);
            const svuint64_t sum = svadd_u64_x(pg_all, guarded_x, y);
            const svuint64_t difference = svsub_u64_x(pg_all, svadd_u64_x(pg_all, guarded_x, twice_modulus_vec), y);

            // Apply the inverse-degree scaling to the upper branch.
            const svuint64_t result_x_lazy = MultiplyUIntModLazySVE(pg_all, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);

             // Apply the scaled inverse root to the lower branch.
            const svuint64_t result_y_lazy = MultiplyUIntModLazySVE(pg_all, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            // Normalize the lazy results to [0, q).
            const svuint64_t result_x = ReduceModFactor2To1SVE(pg_all, result_x_lazy, modulus_vec);
            const svuint64_t result_y = ReduceModFactor2To1SVE(pg_all, result_y_lazy, modulus_vec);

            svst1_u64(pg_all, x_ptr + j, result_x);
            svst1_u64(pg_all, y_ptr + j, result_y);
        }

         // Handle the remaining elements with predication.
        if (j < gap) {
            const svbool_t pg_tail = svwhilelt_b64(static_cast<uint64_t>(j), static_cast<uint64_t>(gap));
            const svuint64_t x = svld1_u64(pg_tail, x_ptr + j);
            const svuint64_t y = svld1_u64(pg_tail, y_ptr + j);

            // Guard x to the [0, 2q) range.
            const svuint64_t guarded_x = ReduceModFactor4To2SVE(pg_tail, x, twice_modulus_vec);
            const svuint64_t sum = svadd_u64_x(pg_tail, guarded_x, y);
            const svuint64_t difference = svsub_u64_x(pg_tail, svadd_u64_x(pg_tail, guarded_x, twice_modulus_vec), y);

            // Apply the inverse-degree scaling to the upper branch.
            const svuint64_t result_x_lazy = MultiplyUIntModLazySVE(pg_tail, sum, inv_degree_vec, inv_degree_quotient_vec, modulus_vec);

            // Apply the scaled inverse root to the lower branch.
            const svuint64_t result_y_lazy = MultiplyUIntModLazySVE(pg_tail, difference, scaled_root_vec, scaled_root_quotient_vec, modulus_vec);

            // Normalize the lazy results to [0, q).
            const svuint64_t result_x = ReduceModFactor2To1SVE(pg_tail, result_x_lazy, modulus_vec);
            const svuint64_t result_y = ReduceModFactor2To1SVE(pg_tail, result_y_lazy, modulus_vec);

            svst1_u64(pg_tail, x_ptr + j, result_x);
            svst1_u64(pg_tail, y_ptr + j, result_y);
        }
    }


    /// @brief Computes the inverse NTT using SVE with a compile-time or runtime vector length.
    /// @note Lanes = 0 selects the runtime SVE vector length.
    template <size_t Lanes>
    inline void InverseTransformSVE(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo) noexcept
    {
        
        const size_t lanes = Lanes != 0 ? Lanes : static_cast<size_t>(svcntd());

        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        /*
        * Process all stages except the final inverse stage.
        * Generic inverse path:
        *
        *   gap = 1, 2, 4, ...
        *   m   = N/2, N/4, N/8, ...
        */
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
 
        // Root used by the final inverse stage.
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        // Combine the final inverse root with the inverse degree modulo q.
        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{ scaled_root_operand, ComputeShoupQuotient(scaled_root_operand, modulus) };

        // Process the final inverse stage.
        InverseFinalStageSVE(operand, gap, modulus, twice_modulus, scaled_root, inv_degree_modulo);
    }


    /// @brief Computes the inverse NTT specialized for SVE VL=128.
    /// @note Uses two 64-bit SVE lanes and fuses the initial gap = 1 and gap = 2 stages.
    template <>
    inline void InverseTransformSVE<2>(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo) noexcept
    {
        
        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        if (coeff_count >= 8) {
            // Process the first two stages with the fused VL=128 kernel.
            InverseInitialTwoStagesSVE128(operand, coeff_count, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 2;
            gap <<= 2;
        }
        else if (coeff_count >= 4) {
            // Process the initial gap = 1 stage.
            InverseInitialStageSVE128(operand, m, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 1;
            gap <<= 1;
        }

        // Process the remaining inverse stages.
        while (m > 1) {
            InverseStageSVE(operand, m, gap, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 1;
            gap <<= 1;
        }

        /*
        * Same final preparation used by the generic path.
        */
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{ scaled_root_operand, ComputeShoupQuotient(scaled_root_operand, modulus) };

        InverseFinalStageSVE(operand, gap, modulus, twice_modulus, scaled_root, inv_degree_modulo);
    }


    /// @brief Computes the normalized inverse NTT using SVE with a compile-time or runtime vector length.
    /// @note Lanes = 0 selects the runtime SVE vector length.
    template <size_t Lanes>
    inline void InverseTransformNormalizedSVE(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo) noexcept
    {
        const size_t lanes = Lanes != 0 ? Lanes : static_cast<size_t>(svcntd());

        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        // Process all stages except the final inverse stage.
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

        // Root used by the final inverse stage.
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

         // Combine the final inverse root with the inverse degree modulo q.
        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{ scaled_root_operand, ComputeShoupQuotient(scaled_root_operand, modulus) };

        // Process and normalize the final inverse stage.
        InverseFinalStageNormalizedSVE(operand, gap, modulus, twice_modulus, scaled_root, inv_degree_modulo);
    }


    /// @brief Computes the normalized inverse NTT specialized for SVE VL=128.
    /// @note Uses two 64-bit SVE lanes and fuses the initial gap = 1 and gap = 2 stages.
    template <>
    inline void InverseTransformNormalizedSVE<2>(
        uint64_t *__restrict operand,
        size_t coeff_count,
        uint64_t modulus,
        uint64_t twice_modulus,
        const NTTMultiplyOperand *__restrict inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo) noexcept
    {
        size_t root_index = 1;
        size_t m = coeff_count >> 1;
        size_t gap = 1;

        if (coeff_count >= 8) {
            // Process the first two stages with the fused VL=128 kernel.
            InverseInitialTwoStagesSVE128(operand, coeff_count, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 2;
            gap <<= 2;
        }
        else if (coeff_count >= 4) {
            // Process the initial gap = 1 stage.
            InverseInitialStageSVE128(operand, m, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 1;
            gap <<= 1;
        }

        // Process the remaining inverse stages.
        while (m > 1) {
            InverseStageSVE(operand, m, gap, modulus, twice_modulus, inv_root_powers, root_index);

            m >>= 1;
            gap <<= 1;
        }

        // Root used by the final inverse stage.
        const NTTMultiplyOperand &final_inv_root = inv_root_powers[root_index];

        // Combine the final inverse root with the inverse degree modulo q.
        uint64_t scaled_root_operand = MultiplyUIntModLazy(final_inv_root.operand, inv_degree_modulo.operand, inv_degree_modulo.quotient, modulus);

        scaled_root_operand = ReduceModFactor2To1Native(scaled_root_operand, modulus);

        const NTTMultiplyOperand scaled_root{ scaled_root_operand, ComputeShoupQuotient(scaled_root_operand, modulus) };

        // Process and normalize the final inverse stage.
        InverseFinalStageNormalizedSVE(operand, gap, modulus, twice_modulus, scaled_root, inv_degree_modulo);
    }


} // namespace


    /// @brief Computes the in-place forward negacyclic Harvey NTT in lazy form using SVE.
    /// @post Each output coefficient is in the range [0, 4 * modulus).
    void NTTNegacyclicHarveyLazySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *root_powers)
    {
        
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        const size_t lanes = static_cast<size_t>(svcntd());

        // Dispatch according to the runtime SVE vector length.
        switch (lanes) {
        case 2:
            // SVE VL=128.
            ForwardTransformSVE<2>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;

        default:
            // Vector-length-agnostic fallback for all other SVE vector lengths.
            ForwardTransformSVE<0>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;
        }
    }


    /// @brief Computes the in-place normalized forward negacyclic Harvey NTT using SVE.
    /// @post Each output coefficient is in the range [0, modulus).
    void NTTNegacyclicHarveySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *root_powers)
    {
        NTTNegacyclicHarveyLazySVE(operand, coeff_count_power, modulus, root_powers);

        const size_t coeff_count = size_t{ 1 } << coeff_count_power;

        // Forward lazy output: [0, 4q) -> [0, q).
        EltwiseReduceMod(operand, operand, static_cast<uint64_t>(coeff_count), modulus, 4, 1);
    }

/*
    /// @brief Computes the in-place normalized forward negacyclic Harvey NTT using SVE.
    /// @post Each output coefficient is in the range [0, modulus).
    void NTTNegacyclicHarveySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *root_powers)
    {
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        const size_t lanes = static_cast<size_t>(svcntd());

        // Dispatch according to the runtime SVE vector length.
        switch (lanes) {
        case 2:
            // SVE VL=128.
            ForwardTransformNormalizedSVE<2>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;

        default:
            // Vector-length-agnostic fallback for all other SVE vector lengths.
            ForwardTransformNormalizedSVE<0>(operand, coeff_count, modulus, twice_modulus, root_powers);
            break;
        }
    }
 */   
     
    /// @brief Computes the in-place inverse negacyclic Harvey NTT in lazy form using SVE.
    /// @post Each output coefficient is in the range [0, 2 * modulus).
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

        // Dispatch according to the runtime SVE vector length./
        switch (lanes) {
        case 2:
            // SVE VL=128.
            InverseTransformSVE<2>(operand, coeff_count, modulus, twice_modulus, inv_root_powers, inv_degree_modulo);
            break;
        default:
            // Vector-length-agnostic fallback for all other SVE vector lengths.
            InverseTransformSVE<0>(operand, coeff_count, modulus, twice_modulus, inv_root_powers, inv_degree_modulo);
            break;
        }
    }


    /// @brief Computes the in-place normalized inverse negacyclic Harvey NTT using SVE.
    /// @post Each output coefficient is in the range [0, modulus).
    void InverseNTTNegacyclicHarveySVE(
        uint64_t *operand,
        size_t coeff_count_power,
        uint64_t modulus,
        const NTTMultiplyOperand *inv_root_powers,
        NTTMultiplyOperand inv_degree_modulo)
    {
        const size_t coeff_count = size_t{ 1 } << coeff_count_power;
        const uint64_t twice_modulus = modulus << 1;
        const size_t lanes = static_cast<size_t>(svcntd());

        // Dispatch according to the runtime SVE vector length.
        switch (lanes) {
        case 2:
            // SVE VL=128.
            InverseTransformNormalizedSVE<2>(operand, coeff_count, modulus, twice_modulus, inv_root_powers, inv_degree_modulo);
            break;

        default:
            // Vector-length-agnostic fallback for all other SVE vector lengths.
            InverseTransformNormalizedSVE<0>(operand, coeff_count, modulus, twice_modulus, inv_root_powers, inv_degree_modulo);
            break;
        }
    }


/*
* Test-only SVE hooks.
*/
#if defined(AXHEL_TESTING) && defined(AXHEL_HAS_SVE)

namespace test_hooks {

    // Exposes the SVE forward Harvey butterfly for testing.
    void ForwardButterflySVE(
        svbool_t pg,
        svuint64_t x,
        svuint64_t y,
        svuint64_t root,
        svuint64_t root_quotient,
        svuint64_t modulus,
        svuint64_t twice_modulus,
        svuint64_t& result_x,
        svuint64_t& result_y) noexcept
    {
        ForwardButterflyVector(pg, x, y, root, root_quotient, modulus, twice_modulus, result_x, result_y);
    }


    // Exposes the SVE inverse Harvey butterfly for testing.
    void InverseButterflySVE(
        svbool_t pg,
        svuint64_t x,
        svuint64_t y,
        svuint64_t inv_root,
        svuint64_t inv_root_quotient,
        svuint64_t modulus,
        svuint64_t twice_modulus,
        svuint64_t& result_x,
        svuint64_t& result_y) noexcept
    {
        InverseButterflyVector(pg, x, y, inv_root, inv_root_quotient, modulus, twice_modulus, result_x, result_y);
    }

} // namespace test_hooks
#endif



} // namespace axhel
} // namespace unipi

#endif