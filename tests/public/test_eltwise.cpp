#include <gtest/gtest.h>
#include <axhel/axhel.hpp>

#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace {
using namespace axhel_test;

void expect_eq_vec(const std::vector<std::uint64_t>& got,
                   const std::vector<std::uint64_t>& exp,
                   const char* kernel, std::uint64_t q, std::uint64_t factor) {
    ASSERT_EQ(got.size(), exp.size());
    for (std::size_t i = 0; i < got.size(); ++i) {
        ASSERT_EQ(got[i], exp[i])
            << kernel << " i=" << i << " q=" << q << " factor=" << factor;
    }
}

class GuardedBuffer {
public:
    static constexpr std::size_t guard_words = 8;
    static constexpr std::uint64_t left_canary  = 0x13579bdf2468ace0ULL;
    static constexpr std::uint64_t right_canary = 0xfdb97531eca86420ULL;

    explicit GuardedBuffer(std::size_t n)
        : n_(n), storage_(n + 2 * guard_words, 0) {
        reset_guards();
    }

    std::uint64_t* data() { return storage_.data() + guard_words; }
    const std::uint64_t* data() const { return storage_.data() + guard_words; }
    std::size_t size() const { return n_; }

    void load(const std::vector<std::uint64_t>& values) {
        ASSERT_EQ(values.size(), n_);
        std::copy(values.begin(), values.end(), data());
    }

    std::vector<std::uint64_t> payload() const {
        return std::vector<std::uint64_t>(data(), data() + n_);
    }

    void expect_guards_intact(const char* name, std::size_t n) const {
        for (std::size_t i = 0; i < guard_words; ++i) {
            ASSERT_EQ(storage_[i], left_canary)
                << name << " left guard overwritten, n=" << n << " guard_i=" << i;
            ASSERT_EQ(storage_[guard_words + n_ + i], right_canary)
                << name << " right guard overwritten, n=" << n << " guard_i=" << i;
        }
    }

private:
    void reset_guards() {
        std::fill(storage_.begin(), storage_.begin() + guard_words, left_canary);
        std::fill(storage_.begin() + guard_words + n_, storage_.end(), right_canary);
    }

    std::size_t n_;
    std::vector<std::uint64_t> storage_;
};

TEST(EltwiseAddMod, BoundaryVectorVector) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        const auto vals = boundary_values(q, 1);
        for (auto n : element_sizes()) {
            auto a = make_pattern(n, vals);
            auto b = make_pattern(n, std::vector<std::uint64_t>(vals.rbegin(), vals.rend()));
            std::vector<std::uint64_t> got(n), exp(n);
            for (std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],b[i],q);
            unipi::axhel::EltwiseAddMod(got.data(),a.data(),b.data(),n,q);
            expect_eq_vec(got,exp,"add-vv",q,1);
        }
    }
}

TEST(EltwiseAddMod, VectorScalarCornerCases) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        const auto vals = boundary_values(q, 1);
        for (auto scalar : vals) {
            const std::size_t n = vector_lanes_u64()*2 + 1;
            auto a = make_pattern(n, vals);
            std::vector<std::uint64_t> got(n), exp(n);
            for (std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],scalar,q);
            unipi::axhel::EltwiseAddMod(got.data(),a.data(),scalar,n,q);
            expect_eq_vec(got,exp,"add-vs",q,1);
        }
    }
}

TEST(EltwiseSubMod, BorrowAndNoBorrowInSameVector) {
    for (auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);
        const std::vector<std::uint64_t> va{0,q-1,0,q-1,1,q-2};
        const std::vector<std::uint64_t> vb{q-1,0,1,q-2,2,q-1};
        for (auto n : element_sizes()) {
            auto a=make_pattern(n,va), b=make_pattern(n,vb);
            std::vector<std::uint64_t> got(n), exp(n);
            for(std::size_t i=0;i<n;++i) exp[i]=sub_mod(a[i],b[i],q);
            unipi::axhel::EltwiseSubMod(got.data(),a.data(),b.data(),n,q);
            expect_eq_vec(got,exp,"sub-vv",q,1);
        }
    }
}

TEST(EltwiseSubMod, VectorScalarBoundary) {
    for (auto bits : modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const auto vals=boundary_values(q,1);
        for(auto scalar: vals) {
            auto a=make_pattern(vector_lanes_u64()*2+1,vals);
            std::vector<std::uint64_t> got(a.size()),exp(a.size());
            for(std::size_t i=0;i<a.size();++i) exp[i]=sub_mod(a[i],scalar,q);
            unipi::axhel::EltwiseSubMod(got.data(),a.data(),scalar,a.size(),q);
            expect_eq_vec(got,exp,"sub-vs",q,1);
        }
    }
}


TEST(EltwiseAddMod, ExactCanonicalReductionBoundaries) {
    // These cases exercise the exact transition used by a candidate+umin
    // implementation:
    //   sum < q, sum == q, and q < sum <= 2q-2.
    for (unsigned bits : {20U, 40U, 60U, 62U}) {
        const auto q = find_ntt_prime(1024, bits);
        const std::vector<std::uint64_t> va{
            0, 0, 1, q-1, q-1, q-2, q-2, 1
        };
        const std::vector<std::uint64_t> vb{
            0, q-1, q-1, 1, q-1, 1, 2, 1
        };
        const std::size_t n = 4 * vector_lanes_u64() + 3;
        const auto a = make_pattern(n, va);
        const auto b = make_pattern(n, vb);

        std::vector<std::uint64_t> got(n), exp(n);
        for (std::size_t i = 0; i < n; ++i) exp[i] = add_mod(a[i], b[i], q);
        unipi::axhel::EltwiseAddMod(got.data(), a.data(), b.data(), n, q);
        expect_eq_vec(got, exp, "add exact vv", q, 1);

        const auto scalar_inputs = make_pattern(
            n, std::vector<std::uint64_t>{0,1,q-2,q-1});
        for (std::uint64_t scalar : std::vector<std::uint64_t>{0,1,q-2,q-1}) {
            for (std::size_t i = 0; i < n; ++i)
                exp[i] = add_mod(scalar_inputs[i], scalar, q);
            unipi::axhel::EltwiseAddMod(
                got.data(), scalar_inputs.data(), scalar, n, q);
            expect_eq_vec(got, exp, "add exact vs", q, 1);
        }
    }
}

TEST(EltwiseSubMod, ExactCanonicalWraparoundBoundaries) {
    // These cases directly exercise the uint64_t wraparound used by a
    // candidate+umin implementation.  In particular:
    //   0 - 1       -> UINT64_MAX before correction
    //   0 - (q - 1) -> 2^64 - (q - 1) before correction
    // together with equal and no-borrow cases.
    for (unsigned bits : {20U, 40U, 60U, 62U}) {
        const auto q = find_ntt_prime(1024, bits);
        const std::vector<std::uint64_t> va{
            0, q-1, 0, q-1, 1, q-2, 1, q-1
        };
        const std::vector<std::uint64_t> vb{
            0, q-1, 1, 0, 0, q-1, q-1, 1
        };
        const std::size_t n = 4 * vector_lanes_u64() + 3;
        const auto a = make_pattern(n, va);
        const auto b = make_pattern(n, vb);

        std::vector<std::uint64_t> got(n), exp(n);
        for (std::size_t i = 0; i < n; ++i) exp[i] = sub_mod(a[i], b[i], q);
        unipi::axhel::EltwiseSubMod(got.data(), a.data(), b.data(), n, q);
        expect_eq_vec(got, exp, "sub exact vv", q, 1);

        const auto scalar_inputs = make_pattern(
            n, std::vector<std::uint64_t>{0,1,q-2,q-1});
        for (std::uint64_t scalar : std::vector<std::uint64_t>{0,1,q-2,q-1}) {
            for (std::size_t i = 0; i < n; ++i)
                exp[i] = sub_mod(scalar_inputs[i], scalar, q);
            unipi::axhel::EltwiseSubMod(
                got.data(), scalar_inputs.data(), scalar, n, q);
            expect_eq_vec(got, exp, "sub exact vs", q, 1);
        }
    }
}

TEST(EltwiseAddSubMod, ExhaustiveSmallModuli) {
    // Exhaustively verify every canonical (a,b) pair for small moduli.
    // The reference functions use an independent wide-arithmetic oracle.
    for (std::uint64_t q : {2ULL, 3ULL, 5ULL, 17ULL, 257ULL}) {
        const std::size_t n = static_cast<std::size_t>(q * q);
        std::vector<std::uint64_t> a(n), b(n);
        std::size_t k = 0;
        for (std::uint64_t av = 0; av < q; ++av) {
            for (std::uint64_t bv = 0; bv < q; ++bv, ++k) {
                a[k] = av;
                b[k] = bv;
            }
        }

        std::vector<std::uint64_t> add_got(n), sub_got(n);
        unipi::axhel::EltwiseAddMod(add_got.data(), a.data(), b.data(), n, q);
        unipi::axhel::EltwiseSubMod(sub_got.data(), a.data(), b.data(), n, q);

        for (std::size_t i = 0; i < n; ++i) {
            ASSERT_EQ(add_got[i], add_mod(a[i], b[i], q))
                << "exhaustive add-vv q=" << q << " i=" << i
                << " a=" << a[i] << " b=" << b[i];
            ASSERT_EQ(sub_got[i], sub_mod(a[i], b[i], q))
                << "exhaustive sub-vv q=" << q << " i=" << i
                << " a=" << a[i] << " b=" << b[i];
            ASSERT_LT(add_got[i], q);
            ASSERT_LT(sub_got[i], q);
        }

        // Exercise the separate vector-scalar overloads exhaustively too.
        std::vector<std::uint64_t> all_a(static_cast<std::size_t>(q));
        for (std::uint64_t av = 0; av < q; ++av)
            all_a[static_cast<std::size_t>(av)] = av;

        std::vector<std::uint64_t> add_vs(all_a.size()), sub_vs(all_a.size());
        for (std::uint64_t scalar = 0; scalar < q; ++scalar) {
            unipi::axhel::EltwiseAddMod(
                add_vs.data(), all_a.data(), scalar, all_a.size(), q);
            unipi::axhel::EltwiseSubMod(
                sub_vs.data(), all_a.data(), scalar, all_a.size(), q);

            for (std::size_t i = 0; i < all_a.size(); ++i) {
                ASSERT_EQ(add_vs[i], add_mod(all_a[i], scalar, q))
                    << "exhaustive add-vs q=" << q
                    << " scalar=" << scalar << " i=" << i;
                ASSERT_EQ(sub_vs[i], sub_mod(all_a[i], scalar, q))
                    << "exhaustive sub-vs q=" << q
                    << " scalar=" << scalar << " i=" << i;
                ASSERT_LT(add_vs[i], q);
                ASSERT_LT(sub_vs[i], q);
            }
        }
    }
}

TEST(EltwiseAddSubMod, AllLengthsThroughMaxSVEILP4Boundary) {
    // SVE supports up to 2048-bit vectors => at most 32 uint64_t lanes.
    // Testing every n in [0,129] therefore crosses L and 4L boundaries for
    // every legal SVE vector length, including the 4*32=128 ILP4 boundary.
    const auto q = find_ntt_prime(1024, 62);
    const std::vector<std::uint64_t> va{0,1,q-2,q-1,0,q-1,1,q-2};
    const std::vector<std::uint64_t> vb{q-1,0,1,q-2,1,q-1,q-1,2};

    for (std::size_t n = 0; n <= 129; ++n) {
        const auto a = make_pattern(n, va);
        const auto b = make_pattern(n, vb);
        std::vector<std::uint64_t> got(n), exp(n);

        for (std::size_t i = 0; i < n; ++i) exp[i] = add_mod(a[i], b[i], q);
        unipi::axhel::EltwiseAddMod(got.data(), a.data(), b.data(), n, q);
        expect_eq_vec(got, exp, "add all lengths vv", q, 1);

        for (std::size_t i = 0; i < n; ++i) exp[i] = sub_mod(a[i], b[i], q);
        unipi::axhel::EltwiseSubMod(got.data(), a.data(), b.data(), n, q);
        expect_eq_vec(got, exp, "sub all lengths vv", q, 1);

        // q-1 is a useful scalar for both kernels: ADD crosses the reduction
        // threshold for every nonzero a, while SUB mixes borrow/no-borrow.
        const std::uint64_t scalar = q - 1;

        for (std::size_t i = 0; i < n; ++i) exp[i] = add_mod(a[i], scalar, q);
        unipi::axhel::EltwiseAddMod(got.data(), a.data(), scalar, n, q);
        expect_eq_vec(got, exp, "add all lengths vs", q, 1);

        for (std::size_t i = 0; i < n; ++i) exp[i] = sub_mod(a[i], scalar, q);
        unipi::axhel::EltwiseSubMod(got.data(), a.data(), scalar, n, q);
        expect_eq_vec(got, exp, "sub all lengths vs", q, 1);
    }
}

TEST(EltwiseAddSubMod, WidestSupportedBitWidthBoundaries) {
    // The regression matrix and AXHEL benchmark currently exercise moduli
    // through 62 bits.  Keep a dedicated test at that widest supported width
    // so future arithmetic rewrites cannot silently break the high-end case.
    const auto q = find_ntt_prime(1024, 62);
    const std::size_t n = 4 * vector_lanes_u64() + 1;

    const auto a = make_pattern(
        n, std::vector<std::uint64_t>{0,1,q-2,q-1,q-1,0});
    const auto b = make_pattern(
        n, std::vector<std::uint64_t>{q-1,q-1,2,1,q-1,1});

    std::vector<std::uint64_t> add_got(n), sub_got(n);
    unipi::axhel::EltwiseAddMod(add_got.data(), a.data(), b.data(), n, q);
    unipi::axhel::EltwiseSubMod(sub_got.data(), a.data(), b.data(), n, q);

    for (std::size_t i = 0; i < n; ++i) {
        ASSERT_EQ(add_got[i], add_mod(a[i], b[i], q)) << "i=" << i;
        ASSERT_EQ(sub_got[i], sub_mod(a[i], b[i], q)) << "i=" << i;
        ASSERT_LT(add_got[i], q);
        ASSERT_LT(sub_got[i], q);
    }
}

TEST(EltwiseAddSubMod, DeterministicRandom) {
    std::mt19937_64 rng(test_seed() ^ 0x414444535542ULL);
    const auto sizes=element_sizes();
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        for(int r=0;r<random_repetitions();++r) {
            const auto n=sizes[static_cast<std::size_t>(r)%sizes.size()];
            auto a=random_vector(n,q,rng),b=random_vector(n,q,rng);
            std::vector<std::uint64_t> add_got(n),sub_got(n);
            unipi::axhel::EltwiseAddMod(add_got.data(),a.data(),b.data(),n,q);
            unipi::axhel::EltwiseSubMod(sub_got.data(),a.data(),b.data(),n,q);
            for(std::size_t i=0;i<n;++i) {
                ASSERT_EQ(add_got[i],add_mod(a[i],b[i],q));
                ASSERT_EQ(sub_got[i],sub_mod(a[i],b[i],q));
            }
        }
    }
}

TEST(EltwiseMulMod, BoundaryAllModFactorsAndBarrettBitWidths) {
    for(const auto& mc: barrett_modulus_cases()) {
        const auto q=mc.q;
        for(std::uint64_t factor: {1ULL,2ULL,4ULL}) {
            const auto vals=boundary_values(q,factor);
            for(auto n: element_sizes()) {
                auto a=make_pattern(n,vals);
                auto b=make_pattern(n,std::vector<std::uint64_t>(vals.rbegin(),vals.rend()));
                std::vector<std::uint64_t> got(n),exp(n);
                for(std::size_t i=0;i<n;++i) exp[i]=mul_mod(a[i]%q,b[i]%q,q);
                unipi::axhel::EltwiseMulMod(got.data(),a.data(),b.data(),n,q,factor);
                expect_eq_vec(got,exp,"mul-vv",q,factor);
            }
        }
    }
}

TEST(EltwiseMulMod, QMinusOneSquared) {
    for(auto bits: modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        const std::size_t n=vector_lanes_u64()*3+1;
        std::vector<std::uint64_t> a(n,q-1),b(n,q-1),got(n);
        unipi::axhel::EltwiseMulMod(got.data(),a.data(),b.data(),n,q,1);
        for(auto x:got) EXPECT_EQ(x,1ULL) << "q="<<q;
    }
}

TEST(EltwiseMulMod, DeterministicRandomAllBarrettBitWidths) {
    std::mt19937_64 rng(test_seed());
    const auto sizes=element_sizes();
    for(const auto& mc: barrett_modulus_cases()) {
        const auto q=mc.q;
        for(std::uint64_t factor:{1ULL,2ULL,4ULL}) {
            const auto bound=checked_bound(q,factor);
            for(int r=0;r<random_repetitions();++r) {
                const auto n=sizes[static_cast<std::size_t>(r)%sizes.size()];
                auto a=random_vector(n,bound,rng),b=random_vector(n,bound,rng);
                std::vector<std::uint64_t> got(n);
                unipi::axhel::EltwiseMulMod(got.data(),a.data(),b.data(),n,q,factor);
                for(std::size_t i=0;i<n;++i)
                    ASSERT_EQ(got[i],mul_mod(a[i]%q,b[i]%q,q))
                        <<"seed="<<test_seed()<<" r="<<r<<" i="<<i
                        <<" bits="<<mc.bits<<" q="<<q<<" q_variant="<<mc.variant
                        <<" factor="<<factor;
            }
        }
    }
}

TEST(EltwiseFMAMod, BoundaryAllModFactorsAndBarrettBitWidths) {
    for(const auto& mc: barrett_modulus_cases()) {
        const auto q=mc.q;
        for(std::uint64_t factor:{1ULL,2ULL,4ULL,8ULL}) {
            const auto vals=boundary_values_uint64_domain(q,factor);
            for(auto n: element_sizes()) {
                auto a=make_pattern(n,vals);
                auto c=make_pattern(n,std::vector<std::uint64_t>(vals.rbegin(),vals.rend()));
                const auto scalar=vals[vals.size()/2];
                std::vector<std::uint64_t> got(n),exp(n);
                for(std::size_t i=0;i<n;++i)
                    exp[i]=fma_mod(a[i]%q,scalar%q,c[i]%q,q);
                unipi::axhel::EltwiseFMAMod(got.data(),a.data(),scalar,c.data(),n,q,factor);
                expect_eq_vec(got,exp,"fma-vsv",q,factor);
            }
        }
    }
}

TEST(EltwiseFMAMod, NullAddendMatchesMultiplication) {
    for(auto bits:modulus_bits()) {
        const auto q=find_ntt_prime(1024,bits);
        for(std::uint64_t factor:{1ULL,2ULL,4ULL,8ULL}) {
            const auto vals=boundary_values_uint64_domain(q,factor);
            const std::size_t n=vector_lanes_u64()*3+1;
            auto a=make_pattern(n,vals);
            const auto scalar=vals.back();
            std::vector<std::uint64_t> got(n);
            unipi::axhel::EltwiseFMAMod(got.data(),a.data(),scalar,nullptr,n,q,factor);
            for(std::size_t i=0;i<n;++i)
                ASSERT_EQ(got[i],mul_mod(a[i]%q,scalar%q,q))
                    <<"bits="<<bits<<" factor="<<factor<<" i="<<i;
        }
    }
}

TEST(EltwiseFMAMod, DeterministicRandomAllBarrettBitWidths) {
    std::mt19937_64 rng(test_seed() ^ 0x464d4152414e44ULL);
    const auto sizes=element_sizes();
    for(const auto& mc: barrett_modulus_cases()) {
        const auto q=mc.q;
        for(std::uint64_t factor:{1ULL,2ULL,4ULL,8ULL}) {
            const auto effective_bound=effective_uint64_bound(q,factor);
            for(int r=0;r<random_repetitions();++r) {
                const auto n=sizes[static_cast<std::size_t>(r)%sizes.size()];
                auto a=random_vector(n,effective_bound,rng);
                auto c=random_vector(n,effective_bound,rng);
                const auto scalar=bounded_random(rng,effective_bound);
                std::vector<std::uint64_t> got(n);
                unipi::axhel::EltwiseFMAMod(got.data(),a.data(),scalar,c.data(),n,q,factor);
                for(std::size_t i=0;i<n;++i)
                    ASSERT_EQ(got[i],fma_mod(a[i]%q,scalar%q,c[i]%q,q))
                        <<"seed="<<test_seed()<<" bits="<<mc.bits<<" q="<<q
                        <<" q_variant="<<mc.variant<<" factor="<<factor
                        <<" r="<<r<<" i="<<i;
            }
        }
    }
}

TEST(EltwiseOperations, InPlaceAliasing) {
    const auto q=find_ntt_prime(1024,50);
    const std::size_t n=vector_lanes_u64()*3+1;
    const auto vals=boundary_values(q,1);
    auto a=make_pattern(n,vals);
    auto b=make_pattern(n,std::vector<std::uint64_t>(vals.rbegin(),vals.rend()));

    // ADD vector-vector: result may alias either input.
    {
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],b[i],q);
        unipi::axhel::EltwiseAddMod(in.data(),in.data(),b.data(),n,q);
        expect_eq_vec(in,exp,"add alias result=op1",q,1);
    }
    {
        auto in=b, exp=b;
        for(std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],b[i],q);
        unipi::axhel::EltwiseAddMod(in.data(),a.data(),in.data(),n,q);
        expect_eq_vec(in,exp,"add alias result=op2",q,1);
    }

    // SUB vector-vector: verify both alias directions independently.
    {
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=sub_mod(a[i],b[i],q);
        unipi::axhel::EltwiseSubMod(in.data(),in.data(),b.data(),n,q);
        expect_eq_vec(in,exp,"sub alias result=op1",q,1);
    }
    {
        auto in=b, exp=b;
        for(std::size_t i=0;i<n;++i) exp[i]=sub_mod(a[i],b[i],q);
        unipi::axhel::EltwiseSubMod(in.data(),a.data(),in.data(),n,q);
        expect_eq_vec(in,exp,"sub alias result=op2",q,1);
    }

    // Same buffer as result and both vector operands.
    {
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],a[i],q);
        unipi::axhel::EltwiseAddMod(in.data(),in.data(),in.data(),n,q);
        expect_eq_vec(in,exp,"add alias all same",q,1);
    }
    {
        auto in=a;
        unipi::axhel::EltwiseSubMod(in.data(),in.data(),in.data(),n,q);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(in[i],0ULL)<<"sub alias all same i="<<i;
    }

    // Vector-scalar overloads may overwrite operand1 in place.
    {
        const std::uint64_t scalar=q-1;
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=add_mod(a[i],scalar,q);
        unipi::axhel::EltwiseAddMod(in.data(),in.data(),scalar,n,q);
        expect_eq_vec(in,exp,"add-vs alias result=op1",q,1);
    }
    {
        const std::uint64_t scalar=q-1;
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=sub_mod(a[i],scalar,q);
        unipi::axhel::EltwiseSubMod(in.data(),in.data(),scalar,n,q);
        expect_eq_vec(in,exp,"sub-vs alias result=op1",q,1);
    }

    {
        auto in=a, exp=a;
        for(std::size_t i=0;i<n;++i) exp[i]=mul_mod(a[i],b[i],q);
        unipi::axhel::EltwiseMulMod(in.data(),in.data(),b.data(),n,q,1);
        expect_eq_vec(in,exp,"mul alias result=op1",q,1);
    }
    {
        auto in=b, exp=b;
        for(std::size_t i=0;i<n;++i) exp[i]=mul_mod(a[i],b[i],q);
        unipi::axhel::EltwiseMulMod(in.data(),a.data(),in.data(),n,q,1);
        expect_eq_vec(in,exp,"mul alias result=op2",q,1);
    }
}

TEST(EltwiseOperations, UnalignedOffsets) {
    const auto q=find_ntt_prime(1024,50);
    const auto L=vector_lanes_u64();
    const std::size_t n=3*L+1;

    for(std::size_t off=0;off<4;++off) {
        std::vector<std::uint64_t> sa(n+8),sb(n+8),so(n+8);
        for(std::size_t i=0;i<n;++i) {
            sa[off+i]=(i*17+3)%q;
            sb[off+i]=(i*29+5)%q;
        }

        unipi::axhel::EltwiseAddMod(
            so.data()+off,sa.data()+off,sb.data()+off,n,q);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(so[off+i],add_mod(sa[off+i],sb[off+i],q))
                <<"add off="<<off<<" i="<<i;

        unipi::axhel::EltwiseSubMod(
            so.data()+off,sa.data()+off,sb.data()+off,n,q);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(so[off+i],sub_mod(sa[off+i],sb[off+i],q))
                <<"sub off="<<off<<" i="<<i;

        const std::uint64_t scalar=q-1;
        unipi::axhel::EltwiseAddMod(
            so.data()+off,sa.data()+off,scalar,n,q);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(so[off+i],add_mod(sa[off+i],scalar,q))
                <<"add-vs off="<<off<<" i="<<i;

        unipi::axhel::EltwiseSubMod(
            so.data()+off,sa.data()+off,scalar,n,q);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(so[off+i],sub_mod(sa[off+i],scalar,q))
                <<"sub-vs off="<<off<<" i="<<i;

        unipi::axhel::EltwiseMulMod(
            so.data()+off,sa.data()+off,sb.data()+off,n,q,1);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(so[off+i],mul_mod(sa[off+i],sb[off+i],q))
                <<"mul off="<<off<<" i="<<i;
    }
}

TEST(EltwiseOperations, GuardCanariesAcrossSVETails) {
    // A 62-bit modulus stresses the largest supported HE Barrett shift while
    // canaries detect any store before/after the requested output interval.
    const auto q=find_ntt_prime(1024,62);

    for(auto n:canary_tail_sizes()) {
        // Add/Sub on canonical inputs.
        const auto vals1=boundary_values(q,1);
        const auto av=make_pattern(n,vals1);
        const auto bv=make_pattern(n,std::vector<std::uint64_t>(vals1.rbegin(),vals1.rend()));
        const auto av_before=av;
        const auto bv_before=bv;

        GuardedBuffer a(n),b(n),out(n);
        a.load(av); b.load(bv);
        unipi::axhel::EltwiseAddMod(out.data(),a.data(),b.data(),n,q);
        a.expect_guards_intact("add op1",n); b.expect_guards_intact("add op2",n);
        out.expect_guards_intact("add result",n);
        EXPECT_EQ(a.payload(),av_before); EXPECT_EQ(b.payload(),bv_before);
        for(std::size_t i=0;i<n;++i) ASSERT_EQ(out.data()[i],add_mod(av[i],bv[i],q));

        GuardedBuffer sub_out(n);
        unipi::axhel::EltwiseSubMod(sub_out.data(),a.data(),b.data(),n,q);
        a.expect_guards_intact("sub op1",n); b.expect_guards_intact("sub op2",n);
        sub_out.expect_guards_intact("sub result",n);
        for(std::size_t i=0;i<n;++i) ASSERT_EQ(sub_out.data()[i],sub_mod(av[i],bv[i],q));

        // Exercise the separate vector-scalar overloads across the same tails.
        const std::uint64_t addsub_scalar=q-1;

        GuardedBuffer add_vs_out(n);
        unipi::axhel::EltwiseAddMod(
            add_vs_out.data(),a.data(),addsub_scalar,n,q);
        a.expect_guards_intact("add-vs op1",n);
        add_vs_out.expect_guards_intact("add-vs result",n);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(add_vs_out.data()[i],add_mod(av[i],addsub_scalar,q));

        GuardedBuffer sub_vs_out(n);
        unipi::axhel::EltwiseSubMod(
            sub_vs_out.data(),a.data(),addsub_scalar,n,q);
        a.expect_guards_intact("sub-vs op1",n);
        sub_vs_out.expect_guards_intact("sub-vs result",n);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(sub_vs_out.data()[i],sub_mod(av[i],addsub_scalar,q));

        // Mul on the widest representable lazy input domain [0,4q).
        const auto vals4=boundary_values(q,4);
        const auto m1=make_pattern(n,vals4);
        const auto m2=make_pattern(n,std::vector<std::uint64_t>(vals4.rbegin(),vals4.rend()));
        GuardedBuffer ma(n),mb(n),mout(n);
        ma.load(m1); mb.load(m2);
        unipi::axhel::EltwiseMulMod(mout.data(),ma.data(),mb.data(),n,q,4);
        ma.expect_guards_intact("mul op1",n); mb.expect_guards_intact("mul op2",n);
        mout.expect_guards_intact("mul result",n);
        EXPECT_EQ(ma.payload(),m1); EXPECT_EQ(mb.payload(),m2);
        for(std::size_t i=0;i<n;++i) ASSERT_EQ(mout.data()[i],mul_mod(m1[i]%q,m2[i]%q,q));

        // FMA factor 8 additionally covers the uint64_t-capped lazy domain.
        const auto vals8=boundary_values_uint64_domain(q,8);
        const auto f1=make_pattern(n,vals8);
        const auto f3=make_pattern(n,std::vector<std::uint64_t>(vals8.rbegin(),vals8.rend()));
        const auto scalar=vals8[vals8.size()/2];
        GuardedBuffer fa(n),fc(n),fout(n);
        fa.load(f1); fc.load(f3);
        unipi::axhel::EltwiseFMAMod(fout.data(),fa.data(),scalar,fc.data(),n,q,8);
        fa.expect_guards_intact("fma op1",n); fc.expect_guards_intact("fma op3",n);
        fout.expect_guards_intact("fma result",n);
        EXPECT_EQ(fa.payload(),f1); EXPECT_EQ(fc.payload(),f3);
        for(std::size_t i=0;i<n;++i)
            ASSERT_EQ(fout.data()[i],fma_mod(f1[i]%q,scalar%q,f3[i]%q,q));
    }
}

} // namespace
