#include <gtest/gtest.h>
#include <axhel/axhel.hpp>
#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

namespace {
using namespace axhel_test;
using unipi::axhel::NTT;

std::vector<std::uint64_t> special_input(std::size_t N, std::uint64_t q, int pattern) {
    std::vector<std::uint64_t> v(N,0);
    switch(pattern) {
        case 0: break;
        case 1: if(N) v[0]=1; break;
        case 2: std::fill(v.begin(),v.end(),1); break;
        case 3: std::fill(v.begin(),v.end(),q-1); break;
        case 4:
            for(std::size_t i=0;i<N;++i) v[i]=(i&1)?q-1:0;
            break;
        case 5:
            if(N) v[N-1]=q-1;
            break;
    }
    return v;
}

TEST(NTT, KnownVectorN8Q17) {
    constexpr std::uint64_t N=8;
    constexpr std::uint64_t q=17;
    constexpr std::uint64_t psi=3;
    const std::vector<std::uint64_t> original{1,2,3,4,5,6,7,8};
    const std::vector<std::uint64_t> expected_forward{5,0,13,8,9,11,5,8};

    NTT ntt(N,q,psi);
    auto x=original;
    ntt.ComputeForward(x.data(),x.data(),1,1);
    EXPECT_EQ(x,expected_forward);

    ntt.ComputeInverse(x.data(),x.data(),1,1);
    EXPECT_EQ(x,original);
}

TEST(NTT, RootSanityAllDegrees) {
    for(auto N:ntt_degrees()) {
        for(auto bits: modulus_bits()) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            EXPECT_EQ(pow_mod(psi,N,q),q-1)<<"N="<<N<<" bits="<<bits;
            EXPECT_EQ(pow_mod(psi,2*N,q),1ULL)<<"N="<<N<<" bits="<<bits;
        }
    }
}

TEST(NTT, ForwardNormalizedAndLazyAreCongruent) {
    for(auto N:ntt_degrees()) {
        for(auto bits:modulus_bits()) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            NTT ntt(N,q,psi);
            for(int pattern=0;pattern<6;++pattern) {
                auto a=special_input(N,q,pattern);
                auto norm=a,lazy=a;
                ntt.ComputeForward(norm.data(),norm.data(),1,1);
                ntt.ComputeForward(lazy.data(),lazy.data(),1,4);
                for(std::size_t i=0;i<N;++i) {
                    ASSERT_LT(norm[i],q);
                    ASSERT_LT(lazy[i],checked_bound(q,4));
                    ASSERT_EQ(norm[i],lazy[i]%q)
                      <<"N="<<N<<" bits="<<bits<<" pattern="<<pattern<<" i="<<i;
                }
            }
        }
    }
}

TEST(NTT, ForwardAcceptsLazyFactor4InputAllDegrees) {
    std::mt19937_64 rng(test_seed() ^ 0x46574434494eULL);
    for(auto N:ntt_degrees()) {
        for(auto bits:modulus_bits()) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            NTT ntt(N,q,psi);

            auto canonical=random_vector(N,q,rng);
            auto lazy=canonical;
            for(std::size_t i=0;i<N;++i) {
                const std::uint64_t k=static_cast<std::uint64_t>(i%4);
                lazy[i]+=k*q;
            }

            ntt.ComputeForward(canonical.data(),canonical.data(),1,1);
            ntt.ComputeForward(lazy.data(),lazy.data(),4,1);
            ASSERT_EQ(lazy,canonical)<<"N="<<N<<" bits="<<bits;
        }
    }
}

TEST(NTT, InverseAcceptsLazyFactor2InputAllDegrees) {
    std::mt19937_64 rng(test_seed() ^ 0x494e5632494eULL);
    for(auto N:ntt_degrees()) {
        for(auto bits:modulus_bits()) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            NTT ntt(N,q,psi);

            auto original=random_vector(N,q,rng);
            auto spectrum=original;
            ntt.ComputeForward(spectrum.data(),spectrum.data(),1,1);
            for(std::size_t i=0;i<N;++i) {
                if(i&1U) spectrum[i]+=q;
            }

            ntt.ComputeInverse(spectrum.data(),spectrum.data(),2,1);
            ASSERT_EQ(spectrum,original)<<"N="<<N<<" bits="<<bits;
        }
    }
}

TEST(NTT, RoundTripNormalizedAndLazyPaths) {
    std::mt19937_64 rng(test_seed()^0x4e5454ULL);
    for(auto N:ntt_degrees()) {
        for(auto bits:modulus_bits()) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            NTT ntt(N,q,psi);

            for(int r=0;r<ntt_random_repetitions();++r) {
                auto original=random_vector(N,q,rng);

                auto x=original;
                ntt.ComputeForward(x.data(),x.data(),1,1);
                ntt.ComputeInverse(x.data(),x.data(),1,1);
                ASSERT_EQ(x,original)
                    <<"N="<<N<<" bits="<<bits<<" r="<<r<<" seed="<<test_seed();

                // AXHEL/SEAL-style lazy composition:
                // forward [0,4q) -> reduce [0,2q) -> inverse factor 2.
                x=original;
                ntt.ComputeForward(x.data(),x.data(),1,4);
                for(auto z:x) ASSERT_LT(z,checked_bound(q,4));
                unipi::axhel::EltwiseReduceMod(x.data(),x.data(),N,q,4,2);
                for(auto z:x) ASSERT_LT(z,checked_bound(q,2));
                ntt.ComputeInverse(x.data(),x.data(),2,1);
                ASSERT_EQ(x,original)
                    <<"lazy roundtrip N="<<N<<" bits="<<bits<<" r="<<r;

                x=original;
                ntt.ComputeForward(x.data(),x.data(),1,1);
                ntt.ComputeInverse(x.data(),x.data(),1,2);
                for(std::size_t i=0;i<N;++i) {
                    ASSERT_LT(x[i],checked_bound(q,2));
                    ASSERT_EQ(x[i]%q,original[i]);
                }
            }
        }
    }
}

TEST(NTT, CriticalModulusVariantsRoundTrip) {
    // Multiple deterministic q values catch q-dependent Shoup/Harvey behavior
    // without multiplying the full NTT matrix by three.
    constexpr std::size_t N=1024;
    std::mt19937_64 rng(test_seed()^0x4e545451564152ULL);
    for(const auto& mc:barrett_modulus_cases()) {
        if(!is_critical_barrett_bits(mc.bits)) continue;
        const auto q=mc.q;
        const auto psi=find_primitive_2n_root(N,q);
        NTT ntt(N,q,psi);
        for(int r=0;r<4;++r) {
            auto original=random_vector(N,q,rng);
            auto x=original;
            ntt.ComputeForward(x.data(),x.data(),1,4);
            unipi::axhel::EltwiseReduceMod(x.data(),x.data(),N,q,4,2);
            ntt.ComputeInverse(x.data(),x.data(),2,1);
            ASSERT_EQ(x,original)
                <<"bits="<<mc.bits<<" q="<<q<<" q_variant="<<mc.variant<<" r="<<r;
        }
    }
}

TEST(NTT, OutOfPlaceMatchesInPlace) {
    for(std::size_t N:{8ULL,128ULL,1024ULL,8192ULL,32768ULL}) {
        const auto q=find_ntt_prime(N,62);
        const auto psi=find_primitive_2n_root(N,q);
        NTT ntt(N,q,psi);
        std::mt19937_64 rng(test_seed()^N);
        auto in=random_vector(N,q,rng);
        auto inplace=in;
        std::vector<std::uint64_t> out(N);

        ntt.ComputeForward(inplace.data(),inplace.data(),1,1);
        ntt.ComputeForward(out.data(),in.data(),1,1);
        EXPECT_EQ(out,inplace)<<"forward N="<<N;

        auto inv_in=inplace;
        ntt.ComputeInverse(inv_in.data(),inv_in.data(),1,1);
        std::vector<std::uint64_t> inv_out(N);
        ntt.ComputeInverse(inv_out.data(),inplace.data(),1,1);
        EXPECT_EQ(inv_out,inv_in)<<"inverse N="<<N;
        EXPECT_EQ(inv_out,in)<<"roundtrip N="<<N;
    }
}

} // namespace
