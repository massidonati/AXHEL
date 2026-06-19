// Copyright (C) University of Pisa - Dept. of Information Engineering
// License-Identifier: Apache-2.0
// Author: M.Donati

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



int main() {

    ExampleEltwiseVectorVectorAddMod();
    ExampleEltwiseVectorScalarAddMod();
    ExampleEltwiseVectorVectorSubMod();
    ExampleEltwiseVectorScalarSubMod();
    ExampleEltwiseMulMod();
    ExampleEltwiseFMAMod();
    ExampleEltwiseReduceMod();

    return 0;
}