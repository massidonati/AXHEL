// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "ntt/ntt-internal.hpp"

#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <cstdint>
#include <random>
#include <vector>

namespace {
using namespace axhel_test;
using unipi::axhel::ComputeShoupQuotient;
using unipi::axhel::MultiplyUIntModLazy;
using unipi::axhel::NTTMultiplyOperand;
using unipi::axhel::ReduceModFactor2To1Native;
using unipi::axhel::ReduceModFactor4To2Native;
using unipi::axhel::ReduceModFactor8To4Native;

TEST(InternalNativeReduction, Factor2To1ExactThresholds) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        for (auto x : std::vector<std::uint64_t>{0,1,q-1,q,q+1,2*q-2,2*q-1}) {
            const auto got = ReduceModFactor2To1Native(x, q);
            ASSERT_LT(got, q);
            ASSERT_EQ(got, x % q) << "q=" << q << " x=" << x;
        }
    }
}

TEST(InternalNativeReduction, Factor8To4ExactThresholds) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        const auto four_q = 4*q;
        const auto vals = boundary_values_uint64_domain(q,8);
        for (auto x : vals) {
            const auto got = ReduceModFactor8To4Native(x, four_q);
            ASSERT_LT(got, four_q);
            ASSERT_EQ(got % q, x % q) << "q=" << q << " x=" << x;
        }
    }
}

TEST(InternalNativeReduction, ReduceInputFactors1_2_4_8) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024,bits);
        for (std::uint64_t factor : {1ULL,2ULL,4ULL,8ULL}) {
            const auto vals = boundary_values_uint64_domain(q,factor);
            for (auto x : vals) {
                std::uint64_t got=0;
                switch(factor) {
                    case 1: got=unipi::axhel::ReduceInputNative<1>(x,q); break;
                    case 2: got=unipi::axhel::ReduceInputNative<2>(x,q); break;
                    case 4: got=unipi::axhel::ReduceInputNative<4>(x,q); break;
                    case 8: got=unipi::axhel::ReduceInputNative<8>(x,q); break;
                }
                ASSERT_LT(got,q);
                ASSERT_EQ(got,x%q)<<"factor="<<factor<<" q="<<q<<" x="<<x;
            }
        }
    }
}

TEST(InternalNativeReduction, Factor4To2ExactThresholds) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        const auto twice_q = 2*q;
        const auto vals=boundary_values(q,4);
        for (auto x : vals) {
            const auto got = ReduceModFactor4To2Native(x, twice_q);
            ASSERT_LT(got, twice_q);
            ASSERT_EQ(got % q, x % q) << "q=" << q << " x=" << x;
        }
    }
}

TEST(InternalShoupNative, QuotientAndLazyMultiplyBoundaries) {
    for (const auto& mc : barrett_modulus_cases()) {
        const auto q = mc.q;
        const std::vector<std::uint64_t> roots{0,1,2,q/2,q-2,q-1};
        const auto xs = boundary_values(q,4);

        for (auto w : roots) {
            const auto wq = ComputeShoupQuotient(w, q);
            for (auto x : xs) {
                const auto got = MultiplyUIntModLazy(x, w, wq, q);
                ASSERT_LT(got, 2*q)
                    << "lazy range bits="<<mc.bits<<" q=" << q << " x=" << x << " w=" << w;
                ASSERT_EQ(got % q, mul_mod(x % q, w, q))
                    << "congruence bits="<<mc.bits<<" q=" << q << " x=" << x << " w=" << w;
            }
        }
    }
}

TEST(InternalShoupNative, DeterministicRandom) {
    std::mt19937_64 rng(test_seed() ^ 0x53484f5550ULL);
    for (const auto& mc : barrett_modulus_cases()) {
        const auto q = mc.q;
        for (int r=0; r<random_repetitions(); ++r) {
            const auto x=bounded_random(rng,4*q);
            const auto w=bounded_random(rng,q);
            const auto wq=ComputeShoupQuotient(w,q);
            const auto got=MultiplyUIntModLazy(x,w,wq,q);
            ASSERT_LT(got,2*q);
            ASSERT_EQ(got%q,mul_mod(x%q,w,q))
                <<"seed="<<test_seed()<<" r="<<r<<" bits="<<mc.bits
                <<" q="<<q<<" q_variant="<<mc.variant;
        }
    }
}

TEST(InternalNTTButterflyNative, ForwardBoundaryAndRange) {
    for (auto bits : modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto twice_q=2*q;
        const auto psi=find_primitive_2n_root(1024,q);
        const NTTMultiplyOperand root{psi,ComputeShoupQuotient(psi,q)};
        const auto vals=boundary_values(q,4);

        for(auto x0:vals) for(auto y0:vals) {
            std::uint64_t x=x0, y=y0;
            unipi::axhel::detail::ForwardButterflyNative(x,y,root,q,twice_q);
            ASSERT_LT(x,4*q);
            ASSERT_LT(y,4*q);
            const auto t=mul_mod(y0%q,psi,q);
            ASSERT_EQ(x%q,add_mod(x0%q,t,q));
            ASSERT_EQ(y%q,sub_mod(x0%q,t,q));
        }
    }
}

TEST(InternalNTTButterflyNative, InverseBoundaryAndRange) {
    for (auto bits : modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto twice_q=2*q;
        const auto psi=find_primitive_2n_root(1024,q);
        const NTTMultiplyOperand inv_root{psi,ComputeShoupQuotient(psi,q)};
        const auto vals=boundary_values(q,2);

        for(auto x0:vals) for(auto y0:vals) {
            std::uint64_t x=x0, y=y0;
            unipi::axhel::detail::InverseButterflyNative(x,y,inv_root,q,twice_q);
            ASSERT_LT(x,2*q);
            ASSERT_LT(y,2*q);
            const auto sum=add_mod(x0%q,y0%q,q);
            const auto diff=sub_mod(x0%q,y0%q,q);
            ASSERT_EQ(x%q,sum);
            ASSERT_EQ(y%q,mul_mod(diff,psi,q));
        }
    }
}

} // namespace
