// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#ifdef AXHEL_HAS_SVE
#include <arm_sve.h>

#include "axhel/number-theory/modular-reduction.hpp"
#include "axhel/number-theory/sve-arith.hpp"
#include "axhel/number-theory/uint-arith.hpp"
#include "ntt/ntt-internal.hpp"

#include "ntt-sve-test-hooks.hpp"
#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#ifndef AXHEL_TEST_LEVEL
#define AXHEL_TEST_LEVEL 0
#endif

namespace {
using namespace axhel_test;
using namespace unipi::axhel;

std::vector<std::uint64_t> store_vec(svbool_t pg, svuint64_t v) {
    std::vector<std::uint64_t> out(svcntd(), 0);
    svst1_u64(pg, out.data(), v);
    return out;
}

TEST(InternalSVEArithmetic, MulU64ToU128MatchesUint128) {
    const auto lanes=static_cast<std::size_t>(svcntd());
    const svbool_t pg=svptrue_b64();

    const std::vector<std::uint64_t> values{
        0,1,2,
        0xffffffffULL,
        0x100000000ULL,
        0x7fffffffffffffffULL,
        0x8000000000000000ULL,
        0xfffffffffffffffeULL,
        0xffffffffffffffffULL
    };

    for(auto a:values) for(auto b:values) {
        const auto va=svdup_n_u64(a), vb=svdup_n_u64(b);
        svuint64_t hi,lo;
        MulU64ToU128SVE(pg,va,vb,&hi,&lo);
        const auto H=store_vec(pg,hi), L=store_vec(pg,lo);
        const u128 p=static_cast<u128>(a)*b;
        const auto exp_lo=static_cast<std::uint64_t>(p);
        const auto exp_hi=static_cast<std::uint64_t>(p>>64);
        for(std::size_t i=0;i<lanes;++i) {
            ASSERT_EQ(L[i],exp_lo)<<"a="<<a<<" b="<<b<<" lane="<<i;
            ASSERT_EQ(H[i],exp_hi)<<"a="<<a<<" b="<<b<<" lane="<<i;
        }
    }
}

TEST(InternalSVEArithmetic, MulHighMatchesUint128) {
    const svbool_t pg=svptrue_b64();
    std::mt19937_64 rng(test_seed()^0x4d554c48494748ULL);
    const int reps=AXHEL_TEST_LEVEL==0?32:(AXHEL_TEST_LEVEL==1?512:4096);
    for(int r=0;r<reps;++r) {
        const auto a=rng(),b=rng();
        const auto got=store_vec(pg,MulHighU64SVE(pg,svdup_n_u64(a),svdup_n_u64(b)));
        const auto exp=static_cast<std::uint64_t>((static_cast<u128>(a)*b)>>64);
        for(auto x:got) ASSERT_EQ(x,exp)<<"r="<<r<<" seed="<<test_seed();
    }
}


TEST(InternalSVEReduction, ConditionalSubtractPreservesInactiveLanes) {
    const auto lanes=static_cast<std::size_t>(svcntd());
    const svbool_t pg_all=svptrue_b64();
    const std::vector<std::uint64_t> bounds{
        1,
        2,
        17,
        std::uint64_t{1} << 62,
        std::numeric_limits<std::uint64_t>::max()
    };

    std::mt19937_64 rng(test_seed()^0x554d494e54455354ULL);

    for(const auto bound:bounds) {
        for(std::size_t active=0;active<=lanes;++active) {
            const svbool_t pg=svwhilelt_b64(
                std::uint64_t{0},static_cast<std::uint64_t>(active));

            for(int repetition=0;repetition<32;++repetition) {
                std::vector<std::uint64_t> input(lanes);
                for(auto& x:input) x=rng();

                if(!input.empty()) input[0]=0;
                if(input.size()>1) input[1]=bound-1;
                if(input.size()>2) input[2]=bound;
                if(input.size()>3) {
                    input[3]=bound==std::numeric_limits<std::uint64_t>::max()
                        ? bound : bound+1;
                }

                const auto x=svld1_u64(pg_all,input.data());
                const auto y=ConditionalSubtractSVE(
                    pg,x,svdup_n_u64(bound));
                const auto output=store_vec(pg_all,y);

                for(std::size_t lane=0;lane<lanes;++lane) {
                    const auto expected=
                        lane<active && input[lane]>=bound
                            ? input[lane]-bound
                            : input[lane];
                    ASSERT_EQ(output[lane],expected)
                        <<"bound="<<bound
                        <<" active="<<active
                        <<" lane="<<lane
                        <<" input="<<input[lane];
                }
            }
        }
    }
}


TEST(InternalSVEReduction, ReduceInputFactor2And4) {
    const svbool_t pg=svptrue_b64();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto vq=svdup_n_u64(q);
        const auto v2q=svdup_n_u64(2*q);
        const auto v4q=svdup_n_u64(4*q);

        for(std::uint64_t factor:{2ULL,4ULL}) {
            const auto vals=boundary_values(q,factor);
            for(std::size_t off=0;off<vals.size();off+=svcntd()) {
                std::vector<std::uint64_t> in(svcntd(),0);
                for(std::size_t i=0;i<in.size();++i)
                    in[i]=vals[(off+i)%vals.size()];
                const auto vin=svld1_u64(pg,in.data());
                svuint64_t vout;
                if(factor==2)
                    vout=ReduceInputSVE<2>(pg,vin,vq,v2q,v4q);
                else
                    vout=ReduceInputSVE<4>(pg,vin,vq,v2q,v4q);
                const auto got=store_vec(pg,vout);
                for(std::size_t i=0;i<got.size();++i) {
                    ASSERT_LT(got[i],q);
                    ASSERT_EQ(got[i],in[i]%q)
                        <<"factor="<<factor<<" bits="<<bits<<" lane="<<i;
                }
            }
        }
    }
}

template <std::uint64_t Shift>
void CheckShiftRight128(std::uint64_t seed, int reps) {
    const svbool_t pg=svptrue_b64();
    std::mt19937_64 rng(seed ^ Shift);
    for(int r=0;r<reps;++r) {
        std::vector<std::uint64_t> hi(svcntd()),lo(svcntd());
        for(auto& x:hi) x=rng();
        for(auto& x:lo) x=rng();
        const auto got=store_vec(
            pg,
            ShiftRight128LowPart<Shift>(
                pg,svld1_u64(pg,hi.data()),svld1_u64(pg,lo.data())));
        for(std::size_t i=0;i<got.size();++i) {
            const u128 x=(static_cast<u128>(hi[i])<<64)|lo[i];
            const auto exp=static_cast<std::uint64_t>(x>>Shift);
            ASSERT_EQ(got[i],exp)
                <<"Shift="<<Shift<<" r="<<r<<" lane="<<i<<" seed="<<seed;
        }
    }
}

TEST(InternalSVEArithmetic, ShiftRight128LowPartRepresentativeBarrettShifts) {
    const int reps=AXHEL_TEST_LEVEL==0?16:(AXHEL_TEST_LEVEL==1?128:1024);
    const auto seed=test_seed()^0x5348494654313238ULL;
    // Representative values cover modulus bit-widths used by the test matrix.
    CheckShiftRight128<18>(seed,reps);
    CheckShiftRight128<28>(seed,reps);
    CheckShiftRight128<38>(seed,reps);
    CheckShiftRight128<48>(seed,reps);
    CheckShiftRight128<50>(seed,reps);
    CheckShiftRight128<53>(seed,reps);
    CheckShiftRight128<57>(seed,reps);
    CheckShiftRight128<58>(seed,reps);
    CheckShiftRight128<59>(seed,reps);
}

TEST(InternalSVEReduction, Factor4To2MixedLanes) {
    const svbool_t pg=svptrue_b64();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto tq=2*q;
        const auto vals=boundary_values(q,4);
        for(std::size_t offset=0;offset<vals.size();offset+=svcntd()) {
            std::vector<std::uint64_t> in(svcntd(),0);
            for(std::size_t i=0;i<in.size();++i)
                in[i]=vals[(offset+i)%vals.size()];
            const auto vin=svld1_u64(pg,in.data());
            const auto got=store_vec(pg,ReduceModFactor4To2SVE(pg,vin,svdup_n_u64(tq)));
            for(std::size_t i=0;i<got.size();++i) {
                const auto expected=in[i]<tq?in[i]:in[i]-tq;
                ASSERT_LT(got[i],tq);
                ASSERT_EQ(got[i]%q,in[i]%q)
                    <<"bits="<<bits<<" lane="<<i<<" input="<<in[i];
                ASSERT_EQ(got[i],expected)
                    <<"exact 4->2 mapping failed bits="<<bits
                    <<" lane="<<i<<" input="<<in[i];
            }
        }
    }
}

TEST(InternalSVEShoup, LazyMultiplyMatchesNativeAndOracle) {
    const svbool_t pg=svptrue_b64();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto vals=boundary_values(q,2);
        for(auto w:std::vector<std::uint64_t>{1,2,q/2,q-2,q-1}) {
            const auto wq=ComputeShoupQuotient(w,q);
            for(std::size_t off=0;off<vals.size();off+=svcntd()) {
                std::vector<std::uint64_t> in(svcntd(),0);
                for(std::size_t i=0;i<in.size();++i)
                    in[i]=vals[(off+i)%vals.size()];
                const auto got=store_vec(
                    pg,
                    MultiplyUIntModLazySVE(
                        pg,
                        svld1_u64(pg,in.data()),
                        svdup_n_u64(w),
                        svdup_n_u64(wq),
                        svdup_n_u64(q)));
                for(std::size_t i=0;i<got.size();++i) {
                    const auto nat=MultiplyUIntModLazy(in[i],w,wq,q);
                    ASSERT_LT(got[i],2*q);
                    ASSERT_EQ(got[i]%q,mul_mod(in[i]%q,w,q));
                    ASSERT_EQ(got[i],nat)
                        <<"Native/SVE representative mismatch bits="<<bits<<" lane="<<i;
                }
            }
        }
    }
}

TEST(InternalSVEButterfly, ForwardMatchesNativeExactly) {
    const svbool_t pg=svptrue_b64();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto tq=2*q;
        const auto w=find_primitive_2n_root(1024,q);
        const auto wq=ComputeShoupQuotient(w,q);
        const auto vals=boundary_values(q,4);

        for(std::size_t off=0;off<vals.size();off+=svcntd()) {
            std::vector<std::uint64_t> xs(svcntd()),ys(svcntd());
            for(std::size_t i=0;i<xs.size();++i) {
                xs[i]=vals[(off+i)%vals.size()];
                ys[i]=vals[(vals.size()-1-((off+i)%vals.size()))];
            }

            svuint64_t rx,ry;
            test_hooks::ForwardButterflySVE(
                pg,svld1_u64(pg,xs.data()),svld1_u64(pg,ys.data()),
                svdup_n_u64(w),svdup_n_u64(wq),
                svdup_n_u64(q),svdup_n_u64(tq),rx,ry);
            const auto X=store_vec(pg,rx),Y=store_vec(pg,ry);

            const NTTMultiplyOperand root{w,wq};
            for(std::size_t i=0;i<X.size();++i) {
                auto nx=xs[i],ny=ys[i];
                detail::ForwardButterflyNative(nx,ny,root,q,tq);
                ASSERT_EQ(X[i],nx)<<"lane="<<i<<" bits="<<bits;
                ASSERT_EQ(Y[i],ny)<<"lane="<<i<<" bits="<<bits;
            }
        }
    }
}

TEST(InternalSVEButterfly, InverseMatchesNativeExactly) {
    const svbool_t pg=svptrue_b64();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto tq=2*q;
        const auto w=find_primitive_2n_root(1024,q);
        const auto wq=ComputeShoupQuotient(w,q);
        const auto vals=boundary_values(q,2);

        for(std::size_t off=0;off<vals.size();off+=svcntd()) {
            std::vector<std::uint64_t> xs(svcntd()),ys(svcntd());
            for(std::size_t i=0;i<xs.size();++i) {
                xs[i]=vals[(off+i)%vals.size()];
                ys[i]=vals[(vals.size()-1-((off+i)%vals.size()))];
            }

            svuint64_t rx,ry;
            test_hooks::InverseButterflySVE(
                pg,svld1_u64(pg,xs.data()),svld1_u64(pg,ys.data()),
                svdup_n_u64(w),svdup_n_u64(wq),
                svdup_n_u64(q),svdup_n_u64(tq),rx,ry);
            const auto X=store_vec(pg,rx),Y=store_vec(pg,ry);

            const NTTMultiplyOperand root{w,wq};
            for(std::size_t i=0;i<X.size();++i) {
                auto nx=xs[i],ny=ys[i];
                detail::InverseButterflyNative(nx,ny,root,q,tq);
                ASSERT_EQ(X[i],nx)<<"lane="<<i<<" bits="<<bits;
                ASSERT_EQ(Y[i],ny)<<"lane="<<i<<" bits="<<bits;
            }
        }
    }
}

} // namespace

#else

TEST(InternalSVE, NotBuiltWithoutSVE) {
    GTEST_SKIP() << "AXHEL_HAS_SVE is not enabled in this build";
}

#endif
