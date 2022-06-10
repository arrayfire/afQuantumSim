/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#include "utils.h"

af::array tensor_product(const af::array& mat1, const af::array& mat2)
{
    af::dim4 dims1 = mat1.dims();
    af::dim4 dims2 = mat2.dims();

    if (mat1.type() != mat2.type() || dims1[2] != 1 || dims1[3] != 1 || dims2[2] != 1 || dims2[3] != 1)
        throw af::exception();

    af::array out = af::tile(mat2, dims1[0], dims1[1]);
    af::array resized_mat1 = af::resize(dims2[0], dims2[1], mat1, AF_INTERP_LOWER);

    return out *= resized_mat1;
}

std::string binary_string(int val, int length)
{
    std::string out;
    out.resize(length + 1);
    
    for (int i = 0; i < length; ++i)
        out[i] = static_cast<bool>(val & (1 << (length - i - 1))) ? '1' : '0';
    
    return out;
}

uint32_t reverse_binary(uint32_t val, int length) noexcept
{
    uint32_t out = val & ~((1 << length) - 1);
    const int half = (length - 1) / 2;
    for (int i = 0; i < length; ++i)
    {
        const int bit = val & (1 << i);
        if (i <= half)
            out |= bit << (length - 1 - i * 2);
        else
            out |= bit >> (i * 2 + 1 - length);

    }
    return out;
}

uint32_t extract_binary(uint32_t val, int first, int last) noexcept
{
    return (val >> first) & ((1 << (1 + last - first)) - 1);
}

std::pair<int64_t, int64_t> approximate_fraction(double f, int64_t md)
{
    std::pair<int64_t, int64_t> result;
	int64_t a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
	int64_t x, d, n = 1;
	int i, neg = 0;
 
	if (md <= 1) 
    { 
        result.second = 1;
        result.first = static_cast<int64_t>(f);
        return result; 
    }
 
	if (f < 0)
    {
        neg = 1;
        f = -f;
    }
 
	while (f != floor(f))
    {
        n <<= 1;
        f *= 2;
    }
	d = f;
 
	for (i = 0; i < 64; i++)
    {
		a = n ? d / n : 0;
		if (i && !a) break;
 
		x = d; d = n; n = x % n;
 
		x = a;
		if (k[1] * a + k[0] >= md) {
			x = (md - k[0]) / k[1];
			if (x * 2 >= a || k[1] >= md)
				i = 65;
			else
				break;
		}
 
		h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
		k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
	}

	result.second = k[1];
	result.first = neg ? -h[1] : h[1];

    return result;
}

int64_t gcd(int64_t a, int64_t b)
{
    for(int64_t temp = a; b != 0; temp = a)
    {
        a = b;
        b = temp % b;
    }
    return a;
}

std::string repeat(int n, const std::string& str)
{
    std::ostringstream oss;
    for(int i = 0; i < n; i++)
        oss << str;
    return oss.str();
}

std::size_t utf8str_len(std::string str)
{
    std::size_t len = 0;
    auto s = str.data();
    while (*s) len += (*s++ & 0xc0) != 0x80;
    return len;
}