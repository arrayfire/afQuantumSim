#pragma once

#include <arrayfire.h>
#include <string>

namespace count{
    template <int... Is>
    struct index {};

    template <int N, int... Is>
    struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

    template <int... Is>
    struct gen_seq<0, Is...> : index<Is...> {};
}

std::string repeat(int n, const std::string& str);
std::size_t utf8str_len(std::string str);

af::array tensor_product(const af::array& arr1, const af::array& arr2);

std::string binary_string(int val, int length);
uint32_t reverse_binary(uint32_t val, int length) noexcept;
uint32_t extract_binary(uint32_t val, int first, int last) noexcept;
std::pair<int64_t, int64_t> approximate_fraction(double value, int64_t max_denominator);
int64_t gcd(int64_t a, int64_t b);

static inline uint32_t fast_pow2(uint32_t pow) { return 1 << pow; }
static inline uint32_t fast_log2(uint32_t val) { int counter = -(val!=0); for(;val;counter++) val >>=1; return counter; }