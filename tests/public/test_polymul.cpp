#include <gtest/gtest.h>
#include <axhel/axhel.hpp>
#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <random>
#include <vector>

namespace {
using namespace axhel_test;
using unipi::axhel::NTT;

void axhel_negacyclic_mul(std::vector<std::uint64_t>& out,
                          const std::vector<std::uint64_t>& a,
                          const std::vector<std::uint64_t>& b,
                          std::uint64_t q, std::uint64_t psi) {
    const auto N=a.size();
    NTT ntt(N,q,psi);
    auto A=a,B=b;
    ntt.ComputeForward(A.data(),A.data(),1,1);
    ntt.ComputeForward(B.data(),B.data(),1,1);
    out.resize(N);
    unipi::axhel::EltwiseMulMod(out.data(),A.data(),B.data(),N,q,1);
    ntt.ComputeInverse(out.data(),out.data(),1,1);
}

TEST(PolynomialMultiply, IdentityZeroAndMonomial) {
    for(std::size_t N:{8ULL,16ULL,32ULL,64ULL}) {
        for(unsigned bits:{40U,62U}) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            std::vector<std::uint64_t> b(N);
            for(std::size_t i=0;i<N;++i) b[i]=(i*17+11)%q;

            std::vector<std::uint64_t> one(N,0),zero(N,0),x(N,0),got;
            one[0]=1; x[1]=1;

            axhel_negacyclic_mul(got,one,b,q,psi);
            EXPECT_EQ(got,b)<<"N="<<N<<" bits="<<bits;

            axhel_negacyclic_mul(got,zero,b,q,psi);
            EXPECT_EQ(got,zero)<<"N="<<N<<" bits="<<bits;

            const auto ref=negacyclic_mul_reference(x,b,q);
            axhel_negacyclic_mul(got,x,b,q,psi);
            EXPECT_EQ(got,ref)<<"N="<<N<<" bits="<<bits;
        }
    }
}

TEST(PolynomialMultiply, IndependentQuadraticOracle) {
    std::mt19937_64 rng(test_seed()^0x504f4c594d554cULL);
    const std::vector<std::size_t> degrees{8,16,32,64};
    const std::vector<unsigned> bits_set{30,40,50,60,62};

    for(auto N:degrees) {
        for(auto bits:bits_set) {
            const auto q=find_ntt_prime(N,bits);
            const auto psi=find_primitive_2n_root(N,q);
            for(int r=0;r<polynomial_oracle_repetitions();++r) {
                auto a=random_vector(N,q,rng);
                auto b=random_vector(N,q,rng);
                const auto ref=negacyclic_mul_reference(a,b,q);
                std::vector<std::uint64_t> got;
                axhel_negacyclic_mul(got,a,b,q,psi);
                ASSERT_EQ(got,ref)
                  <<"N="<<N<<" bits="<<bits<<" r="<<r<<" seed="<<test_seed();
            }
        }
    }
}

TEST(PolynomialMultiply, AllQMinusOne) {
    for(unsigned bits:{60U,62U}) {
        const std::size_t N=32;
        const auto q=find_ntt_prime(N,bits);
        const auto psi=find_primitive_2n_root(N,q);
        std::vector<std::uint64_t> a(N,q-1),b(N,q-1),got;
        const auto ref=negacyclic_mul_reference(a,b,q);
        axhel_negacyclic_mul(got,a,b,q,psi);
        EXPECT_EQ(got,ref)<<"bits="<<bits;
    }
}
} // namespace
