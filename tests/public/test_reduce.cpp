#include <gtest/gtest.h>
#include <axhel/axhel.hpp>
#include "../common/parameters.hpp"
#include "../common/reference.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#ifndef AXHEL_TEST_LEVEL
#define AXHEL_TEST_LEVEL 0
#endif

namespace {
using namespace axhel_test;

void check_reduce_result(
    const std::vector<std::uint64_t>& input,
    const std::vector<std::uint64_t>& output,
    std::uint64_t q,
    std::uint64_t in_factor,
    std::uint64_t out_factor)
{
    ASSERT_EQ(output.size(), input.size());

    const auto range = checked_bound(q, out_factor);
    for (std::size_t i = 0; i < input.size(); ++i) {
        ASSERT_LT(output[i], range)
            << "range invariant failed i=" << i
            << " q=" << q
            << " " << in_factor << "->" << out_factor
            << " input=" << input[i];

        ASSERT_EQ(output[i] % q, input[i] % q)
            << "congruence failed i=" << i
            << " q=" << q
            << " " << in_factor << "->" << out_factor
            << " input=" << input[i];
    }
}


void run_reduce_case(
    std::uint64_t q,
    std::uint64_t in_factor,
    std::uint64_t out_factor,
    std::size_t n,
    const std::vector<std::uint64_t>& values)
{
    const auto input = make_pattern(n, values);

    std::vector<std::uint64_t> output(n);
    unipi::axhel::EltwiseReduceMod(
        output.data(), input.data(), static_cast<std::uint64_t>(n),
        q, in_factor, out_factor);
    check_reduce_result(input, output, q, in_factor, out_factor);

    auto inplace = input;
    unipi::axhel::EltwiseReduceMod(
        inplace.data(), inplace.data(), static_cast<std::uint64_t>(n),
        q, in_factor, out_factor);
    check_reduce_result(input, inplace, q, in_factor, out_factor);
}


void run_identity_case(
    std::uint64_t q,
    std::uint64_t in_factor,
    std::uint64_t out_factor,
    std::size_t n)
{
    const auto input = make_pattern(n, boundary_values(q, in_factor));

    std::vector<std::uint64_t> output(n, std::uint64_t{0xdeadbeef});
    unipi::axhel::EltwiseReduceMod(
        output.data(), input.data(), static_cast<std::uint64_t>(n),
        q, in_factor, out_factor);
    EXPECT_EQ(output, input)
        << "out-of-place identity failed for "
        << in_factor << "->" << out_factor << " n=" << n;

    auto inplace = input;
    unipi::axhel::EltwiseReduceMod(
        inplace.data(), inplace.data(), static_cast<std::uint64_t>(n),
        q, in_factor, out_factor);
    EXPECT_EQ(inplace, input)
        << "in-place identity failed for "
        << in_factor << "->" << out_factor << " n=" << n;
}


std::vector<std::size_t> ilp4_sizes()
{
    const std::size_t lanes = vector_lanes_u64();
    std::vector<std::size_t> sizes{
        0,
        1,
        lanes > 1 ? lanes - 1 : 1,
        lanes,
        lanes + 1,
        4 * lanes - 1,
        4 * lanes,
        4 * lanes + 1,
        8 * lanes - 1,
        8 * lanes,
        8 * lanes + 3
    };

    std::sort(sizes.begin(), sizes.end());
    sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
    return sizes;
}


TEST(EltwiseReduceMod, IdentityExactBothAliases)
{
    const auto q = find_ntt_prime(1024, 50);

    for (const auto& factors : {
             std::pair<std::uint64_t, std::uint64_t>{1, 1},
             {1, 2},
             {2, 2}}) {
        for (const auto n : ilp4_sizes()) {
            run_identity_case(q, factors.first, factors.second, n);
        }
    }
}


TEST(EltwiseReduceMod, SurgicalBoundaries)
{
    for (const auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);

        for (const auto& factors : {
                 std::pair<std::uint64_t, std::uint64_t>{2, 1},
                 {4, 1},
                 {4, 2}}) {
            try {
                static_cast<void>(checked_bound(q, factors.first));
            }
            catch (...) {
                continue;
            }

            const auto values = boundary_values(q, factors.first);
            for (const auto n : element_sizes()) {
                run_reduce_case(
                    q, factors.first, factors.second, n, values);
            }
        }
    }
}


TEST(EltwiseReduceMod, ILP4BoundariesAndTail)
{
    const auto q = find_ntt_prime(1024, 50);

    for (const auto& factors : {
             std::pair<std::uint64_t, std::uint64_t>{2, 1},
             {4, 1},
             {4, 2}}) {
        const auto values = boundary_values(q, factors.first);
        for (const auto n : ilp4_sizes()) {
            run_reduce_case(
                q, factors.first, factors.second, n, values);
        }
    }
}


TEST(EltwiseReduceMod, ExactThresholds)
{
    const auto q = find_ntt_prime(1024, 50);
    const std::vector<std::uint64_t> values{
        0, 1, q - 1, q, q + 1,
        2 * q - 2, 2 * q - 1, 2 * q, 2 * q + 1,
        3 * q - 1, 3 * q, 3 * q + 1,
        4 * q - 2, 4 * q - 1
    };

    run_reduce_case(q, 4, 2, values.size(), values);
    run_reduce_case(q, 4, 1, values.size(), values);
}


#ifdef AXHEL_HAS_SVE
TEST(EltwiseReduceMod, Factor4To2ExactSVERepresentative)
{
    constexpr std::uint64_t q = 17;
    const std::vector<std::uint64_t> input{
        0, 1, q - 1, q, q + 1,
        2 * q - 1, 2 * q, 2 * q + 1,
        3 * q - 1, 3 * q, 3 * q + 1,
        4 * q - 2, 4 * q - 1
    };

    auto expected = input;
    for (auto& x : expected) {
        if (x >= 2 * q) {
            x -= 2 * q;
        }
    }

    std::vector<std::uint64_t> output(input.size());
    unipi::axhel::EltwiseReduceMod(
        output.data(), input.data(),
        static_cast<std::uint64_t>(input.size()), q, 4, 2);
    EXPECT_EQ(output, expected);

    auto inplace = input;
    unipi::axhel::EltwiseReduceMod(
        inplace.data(), inplace.data(),
        static_cast<std::uint64_t>(inplace.size()), q, 4, 2);
    EXPECT_EQ(inplace, expected);
}
#endif


TEST(EltwiseReduceMod, DeterministicRandomAndRangeInvariant)
{
    std::mt19937_64 rng(test_seed() ^ 0x524544554345ULL);

    for (const auto bits : modulus_bits()) {
        const auto q = find_ntt_prime(1024, bits);

        for (const auto& factors : {
                 std::pair<std::uint64_t, std::uint64_t>{2, 1},
                 {4, 1},
                 {4, 2}}) {
            std::uint64_t bound = 0;
            try {
                bound = checked_bound(q, factors.first);
            }
            catch (...) {
                continue;
            }

            const auto sizes = element_sizes();
            for (int repetition = 0;
                 repetition < random_repetitions();
                 ++repetition) {
                const auto n = sizes[
                    static_cast<std::size_t>(repetition) % sizes.size()];
                const auto input = random_vector(n, bound, rng);

                std::vector<std::uint64_t> output(n);
                unipi::axhel::EltwiseReduceMod(
                    output.data(), input.data(),
                    static_cast<std::uint64_t>(n), q,
                    factors.first, factors.second);
                check_reduce_result(
                    input, output, q, factors.first, factors.second);

                auto inplace = input;
                unipi::axhel::EltwiseReduceMod(
                    inplace.data(), inplace.data(),
                    static_cast<std::uint64_t>(n), q,
                    factors.first, factors.second);
                check_reduce_result(
                    input, inplace, q, factors.first, factors.second);
            }
        }
    }
}

} // namespace
