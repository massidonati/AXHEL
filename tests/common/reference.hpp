#pragma once
#include <cstdint>
#include <vector>
#include <limits>
#include <stdexcept>

namespace axhel_test {

using u128 = unsigned __int128;

inline std::uint64_t add_mod(std::uint64_t a, std::uint64_t b, std::uint64_t q) {
    return static_cast<std::uint64_t>((static_cast<u128>(a) + b) % q);
}

inline std::uint64_t sub_mod(std::uint64_t a, std::uint64_t b, std::uint64_t q) {
    return a >= b ? (a - b) % q : static_cast<std::uint64_t>((static_cast<u128>(a) + q - b) % q);
}

inline std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b, std::uint64_t q) {
    return static_cast<std::uint64_t>((static_cast<u128>(a) * b) % q);
}

inline std::uint64_t fma_mod(std::uint64_t a, std::uint64_t b,
                             std::uint64_t c, std::uint64_t q) {
    return static_cast<std::uint64_t>(
        (static_cast<u128>(a) * b + c) % q);
}

inline std::uint64_t pow_mod(std::uint64_t a, std::uint64_t e, std::uint64_t q) {
    std::uint64_t r = 1;
    while (e) {
        if (e & 1U) r = mul_mod(r, a, q);
        a = mul_mod(a, a, q);
        e >>= 1U;
    }
    return r;
}

inline std::vector<std::uint64_t> negacyclic_mul_reference(
    const std::vector<std::uint64_t>& a,
    const std::vector<std::uint64_t>& b,
    std::uint64_t q)
{
    if (a.size() != b.size()) throw std::invalid_argument("size mismatch");
    const std::size_t n = a.size();
    std::vector<std::uint64_t> out(n, 0);

    // O(N^2) oracle. Intended only for small N in regression tests.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const auto p = mul_mod(a[i] % q, b[j] % q, q);
            const std::size_t k = i + j;
            if (k < n) {
                out[k] = add_mod(out[k], p, q);
            } else {
                out[k - n] = sub_mod(out[k - n], p, q); // X^N == -1
            }
        }
    }
    return out;
}

inline void canonicalize(std::vector<std::uint64_t>& v, std::uint64_t q) {
    for (auto& x : v) x %= q;
}

} // namespace axhel_test
