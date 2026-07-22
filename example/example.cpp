// Copyright (C) University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0


#include <cstdint>
#include <iostream>
#include <vector>

#include "axhel/axhel.hpp"

// utility function to check equality of calculated vs expected results
bool CheckEqual(const std::vector<uint64_t>& x,
                const std::vector<uint64_t>& y) {
    if (x.size() != y.size()) {
        std::cout << "Not equal in size\n";
        return false;
    }

    bool is_match = true;

    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] != y[i]) {
            std::cout << "Mismatch at index " << i
                      << " (" << x[i] << " vs " << y[i] << ")\n";
            is_match = false;
        }
    }

    return is_match;
}


// utility function to check congruence of calculated vs expected result in lazy-mode
bool CheckCongruentAndRange(const std::vector<uint64_t>& result, const std::vector<uint64_t>& expected, uint64_t modulus, uint64_t upper_bound) {
    if (result.size() != expected.size()) {
        std::cout << "Not equal in size\n";
        return false;
    }

    bool is_match = true;

    for (size_t i = 0; i < result.size(); ++i) {
        if (result[i] >= upper_bound) {
            std::cout << "Value out of range at index " << i
                      << " (" << result[i]
                      << " >= " << upper_bound << ")\n";
            is_match = false;
        }

        if ((result[i] % modulus) !=
            (expected[i] % modulus)) {
            std::cout << "Mismatch modulo " << modulus
                      << " at index " << i
                      << " (" << result[i]
                      << " vs " << expected[i] << ")\n";
            is_match = false;
        }
    }

    return is_match;
}


void ExampleEltwiseVectorVectorAddMod() {
  std::cout << "Running ExampleEltwiseVectorVectorAddMod...\n";

  std::vector<uint64_t> op1{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<uint64_t> op2{1, 3, 5, 7, 2, 4, 6, 8};
  uint64_t modulus = 10;
  std::vector<uint64_t> exp_out{2, 5, 8, 1, 7, 0, 3, 6};

  unipi::axhel::EltwiseAddMod(op1.data(), op1.data(), op2.data(), op1.size(), modulus);

  CheckEqual(op1, exp_out);
  std::cout << "Done running ExampleEltwiseVectorVectorAddMod\n";
}


void ExampleEltwiseVectorScalarAddMod() {
  std::cout << "Running ExampleEltwiseVectorScalarAddMod...\n";

  std::vector<uint64_t> op1{1, 2, 3, 4, 5, 6, 7, 8};
  uint64_t op2{3};
  uint64_t modulus = 10;
  std::vector<uint64_t> exp_out{4, 5, 6, 7, 8, 9, 0, 1};

  unipi::axhel::EltwiseAddMod(op1.data(), op1.data(), op2, op1.size(), modulus);

  CheckEqual(op1, exp_out);
  std::cout << "Done running ExampleEltwiseVectorScalarAddMod\n";
}


void ExampleEltwiseVectorVectorSubMod() {
  std::cout << "Running ExampleEltwiseVectorVectorSubMod...\n";

  std::vector<uint64_t> op1{7, 4, 1, 0, 9, 5, 3, 2};
  std::vector<uint64_t> op2{2, 6, 1, 9, 4, 5, 8, 3};
  uint64_t modulus = 10;
  std::vector<uint64_t> exp_out{5, 8, 0, 1, 5, 0, 5, 9};

  unipi::axhel::EltwiseSubMod(op1.data(), op1.data(), op2.data(), op1.size(), modulus);

  CheckEqual(op1, exp_out);
  std::cout << "Done running ExampleEltwiseVectorVectorSubMod\n";
}


void ExampleEltwiseVectorScalarSubMod() {
  std::cout << "Running ExampleEltwiseVectorScalarSubMod...\n";

  std::vector<uint64_t> op1{7, 4, 1, 0, 9, 5, 3, 2};
  uint64_t op2{6};
  uint64_t modulus = 10;
  std::vector<uint64_t> exp_out{1, 8, 5, 4, 3, 9, 7, 6};

  unipi::axhel::EltwiseSubMod(op1.data(), op1.data(), op2, op1.size(), modulus);

  CheckEqual(op1, exp_out);
  std::cout << "Done running ExampleEltwiseVectorScalarSubMod\n";
}


void ExampleEltwiseMulMod() {
    std::cout << "Running ExampleEltwiseMulMod...\n";

    {
        std::vector<uint64_t> op1{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<uint64_t> op2{2, 3, 4, 5, 6, 7, 8, 9};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{2, 6, 12, 3, 13, 8, 5, 4};

        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 1;

        unipi::axhel::EltwiseMulMod(res.data(), op1.data(), op2.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    {
        std::vector<uint64_t> op1{16, 15, 14, 13, 12, 11, 10, 9};
        std::vector<uint64_t> op2{16, 15, 14, 13, 12, 11, 10, 9};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{1, 4, 9, 16, 8, 2, 15, 13};

        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 1;

        unipi::axhel::EltwiseMulMod(res.data(), op1.data(), op2.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    {
        std::vector<uint64_t> op1{0, 17, 18, 19, 20, 31, 32, 33};
        std::vector<uint64_t> op2{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{0, 0, 3, 8, 15, 16, 3, 9};

        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 2;

        unipi::axhel::EltwiseMulMod(res.data(), op1.data(), op2.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    {
        std::vector<uint64_t> op1{0, 17, 18, 33, 34, 50, 51, 67};
        std::vector<uint64_t> op2{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{0, 0, 3, 13, 0, 11, 0, 9};

        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 4;

        unipi::axhel::EltwiseMulMod(res.data(), op1.data(), op2.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    std::cout << "Done running ExampleEltwiseMulMod\n";
}


void ExampleEltwiseFMAMod() {
    std::cout << "Running ExampleEltwiseFMAMod...\n";

    {
        std::vector<uint64_t> op1{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<uint64_t> op3{0, 1, 2, 3, 4, 5, 6, 7};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{3, 7, 11, 15, 2, 6, 10, 14};

        const uint64_t op2 = 3;
        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 1;

        unipi::axhel::EltwiseFMAMod(res.data(), op1.data(), op2, op3.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    {
        std::vector<uint64_t> op1{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{3, 6, 9, 12, 15, 1, 4, 7};

        const uint64_t op2 = 3;
        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 1;

        unipi::axhel::EltwiseFMAMod(res.data(), op1.data(), op2, nullptr, op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    {
        std::vector<uint64_t> op1{0, 17, 18, 33, 34, 50, 51, 67};
        std::vector<uint64_t> op3{0, 1, 17, 18, 33, 34, 50, 51};
        std::vector<uint64_t> res(op1.size());

        std::vector<uint64_t> exp_out{0, 1, 3, 15, 16, 14, 16, 14};

        const uint64_t op2 = 3;
        const uint64_t modulus = 17;
        const uint64_t input_mod_factor = 4;

        unipi::axhel::EltwiseFMAMod(res.data(), op1.data(), op2, op3.data(), op1.size(), modulus, input_mod_factor);

        CheckEqual(res, exp_out);
    }

    std::cout << "Done running ExampleEltwiseFMAMod\n";
}


void ExampleEltwiseReduceMod() {
    std::cout << "Running ExampleEltwiseReduceMod...\n";

    {
        std::vector<uint64_t> op{0, 9, 10, 11, 19, 20, 29, 39};
        std::vector<uint64_t> exp_out{0, 9, 0, 1, 9, 0, 9, 9};

        const uint64_t modulus = 10;
        const uint64_t in_mod_factor = 4;
        const uint64_t out_mod_factor = 1;

        unipi::axhel::EltwiseReduceMod(op.data(), op.data(), op.size(), modulus, in_mod_factor, out_mod_factor);

        CheckEqual(op, exp_out);
    }

    {
        std::vector<uint64_t> op{0, 9, 10, 19, 20, 21, 29, 39};
        std::vector<uint64_t> exp_out{0, 9, 10, 19, 10, 11, 19, 19};

        const uint64_t modulus = 10;
        const uint64_t in_mod_factor = 4;
        const uint64_t out_mod_factor = 2;

        unipi::axhel::EltwiseReduceMod(op.data(), op.data(), op.size(), modulus, in_mod_factor, out_mod_factor);

        CheckEqual(op, exp_out);
    }

    {
        std::vector<uint64_t> op{0, 1, 16, 17, 18, 33, 34, 35, 100, 288};
        std::vector<uint64_t> exp_out{0, 1, 16, 0, 1, 16, 0, 1, 15, 16};

        const uint64_t modulus = 17;
        const uint64_t in_mod_factor = modulus;
        const uint64_t out_mod_factor = 1;

        unipi::axhel::EltwiseReduceMod(op.data(), op.data(), op.size(), modulus, in_mod_factor, out_mod_factor);

        CheckEqual(op, exp_out);
    }

    std::cout << "Done running ExampleEltwiseReduceMod\n";
}


void ExampleNTT() {
    std::cout << "Running ExampleNTT...\n";

    constexpr size_t coeff_count_power = 3;
    constexpr uint64_t modulus = 17;
    constexpr uint64_t twice_modulus = 2 * modulus;
    constexpr uint64_t four_times_modulus = 4 * modulus;

    /*
     * Polynomial degree:
     *
     *   N = 2^3 = 8
     */
    const std::vector<uint64_t> original{1, 2, 3, 4, 5, 6, 7, 8};

    /*
     * Forward root powers for N = 8 and q = 17.
     *
     * A primitive 16-th root of unity modulo 17 is 3.
     * The values are stored in the order expected by the
     * Harvey forward NTT.
     */
    const std::vector<unipi::axhel::NTTMultiplyOperand>
        root_powers{
            {1,  1085102592571150095ULL},
            {13, 14106333703424951235ULL},
            {9,  9765923333140350855ULL},
            {15, 16276538888567251425ULL},
            {3,  3255307777713450285ULL},
            {5,  5425512962855750475ULL},
            {10, 10851025925711500950ULL},
            {11, 11936128518282651045ULL}
        };

    /*
     * Inverse root powers in the order expected by the
     * Gentleman-Sande inverse NTT.
     */
    const std::vector<unipi::axhel::NTTMultiplyOperand>
        inv_root_powers{
            {0,  0},
            {6,  6510615555426900570ULL},
            {7,  7595718147998050665ULL},
            {12, 13021231110853801140ULL},
            {14, 15191436295996101330ULL},
            {2,  2170205185142300190ULL},
            {8,  8680820740569200760ULL},
            {4,  4340410370284600380ULL}
        };

    /*
     * 8^-1 mod 17 = 15.
     */
    const unipi::axhel::NTTMultiplyOperand inv_degree_modulo{15, 16276538888567251425ULL};

    /*
     * Expected normalized forward NTT result.
     */
    const std::vector<uint64_t> expected_forward{5, 0, 13, 8, 9, 11, 5, 8};

    /*
     * Normalized forward NTT.
     */
    {
        std::vector<uint64_t> result = original;

        unipi::axhel::NTTNegacyclicHarvey(result.data(), coeff_count_power, modulus, root_powers.data());

        const bool success = CheckEqual(result, expected_forward);

        std::cout << "Forward NTT: " << (success ? "PASS" : "FAIL") << '\n';
    }

    /*
     * Normalized inverse NTT.
     *
     * Start from the normalized forward result and verify that
     * the original polynomial is recovered exactly.
     */
    {
        std::vector<uint64_t> result = expected_forward;

        unipi::axhel::InverseNTTNegacyclicHarvey(result.data(), coeff_count_power, modulus, inv_root_powers.data(), inv_degree_modulo);

        const bool success = CheckEqual(result, original);

        std::cout << "Inverse NTT: " << (success ? "PASS" : "FAIL") << '\n';
    }

    /*
     * Lazy forward NTT.
     *
     * The output is not necessarily equal to the normalized
     * expected result, but each coefficient must:
     *
     *   - be congruent modulo q;
     *   - lie in [0, 4q).
     */
    std::vector<uint64_t> forward_lazy = original;

    unipi::axhel::NTTNegacyclicHarveyLazy( forward_lazy.data(), coeff_count_power, modulus, root_powers.data());

    const bool forward_lazy_success = CheckCongruentAndRange(forward_lazy, expected_forward, modulus, four_times_modulus);

    std::cout << "Forward NTT lazy: " << (forward_lazy_success ? "PASS" : "FAIL") << '\n';

    /*
     * Inverse lazy input must lie in [0, 2q).
     *
     * The forward lazy result lies in [0, 4q), so reduce it
     * from factor 4 to factor 2 before invoking the inverse.
     */
    std::vector<uint64_t> inverse_lazy = forward_lazy;

    for (uint64_t& value : inverse_lazy) {
        if (value >= twice_modulus) {
            value -= twice_modulus;
        }
    }

    /*
     * Lazy inverse NTT.
     *
     * The result must:
     *
     *   - be congruent to the original polynomial modulo q;
     *   - lie in [0, 2q).
     */
    unipi::axhel::InverseNTTNegacyclicHarveyLazy( inverse_lazy.data(), coeff_count_power, modulus, inv_root_powers.data(), inv_degree_modulo);

    const bool inverse_lazy_success = CheckCongruentAndRange( inverse_lazy, original, modulus, twice_modulus);

    std::cout << "Inverse NTT lazy: " << (inverse_lazy_success ? "PASS" : "FAIL") << '\n';

    std::cout << "Done running ExampleNTT\n";
}



int main() {

    ExampleEltwiseVectorVectorAddMod();
    ExampleEltwiseVectorScalarAddMod();
    ExampleEltwiseVectorVectorSubMod();
    ExampleEltwiseVectorScalarSubMod();
    ExampleEltwiseMulMod();
    ExampleEltwiseFMAMod();
    ExampleEltwiseReduceMod();
    ExampleNTT();

    return 0;
}