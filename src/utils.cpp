/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#include "utils.h"

static af::array tensor_product_af_cpu(const af::array& mat1, const af::array& mat2);
static af::array tensor_product_af_opencl(const af::array& mat1, const af::array& mat2);

af::array tensor_product(const af::array& mat1, const af::array& mat2)
{
    if (mat1.type() != mat2.type())
        throw std::invalid_argument{"The tensor requires both matrices to have the same type"};
    if (mat1.dims()[2] != 1 || mat1.dims()[3] != 1 || mat2.dims()[2] != 1|| mat2.dims()[3] != 1)
        throw std::invalid_argument{"Tensor product only supports two dimensions"};

    switch (af::getActiveBackend())
    {
    case AF_BACKEND_CPU:
        return tensor_product_af_cpu(mat1, mat2);
    case AF_BACKEND_OPENCL:
         return tensor_product_af_opencl(mat1, mat2);
    default:
        return tensor_product_af_opencl(mat1, mat2);
    }
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
 
	while (f != std::floor(f))
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

af::array insert_bits(const af::array& current, const af::array& insert, const af::array& position, const af::array& offset, size_t array_length)
{
    if (current.dims()[0] != array_length || 
        insert.dims()[0] != array_length || position.dims()[0] != array_length || offset.dims()[0] != array_length)
        throw af::exception{"Different array lengths"};

    af::array least_mask = ((af::constant(1, array_length, s32) << position) - 1);
    af::array least = current & least_mask;
    af::array temp = insert << position;
    temp = temp | least;

    af::array top_mask = ~least_mask;
    af::array top = (current & top_mask) << offset;
    temp = temp | top;

    return temp;
}

af::array gen_index(uint32_t gate_qubit_begin, uint32_t gate_qubit_count, uint32_t circuit_qubit_count, const af::array& values, const af::array& indices,
                    size_t array_length)
{
    if (values.dims()[0] != array_length || indices.dims()[0] != array_length)
        throw af::exception{"Different array lengths"};
    return insert_bits(indices,
                       values,
                       af::constant(circuit_qubit_count - gate_qubit_count - gate_qubit_begin, array_length, s32),
                       af::constant(gate_qubit_count, array_length, s32),
                       array_length);
}

af::array tensor_product_af_cpu(const af::array& mat1, const af::array& mat2)
{
    const auto rows1 = mat1.dims()[0];
    const auto cols1 = mat1.dims()[1];
    const auto rows2 = mat2.dims()[0];
    const auto cols2 = mat2.dims()[1];

    auto temp = af::moddims(af::flat(af::tile(af::flat(mat1).T(), cols2)).T(), cols1 * cols2, rows1);

    auto resized_mat1 = af::moddims(af::tile(temp.T(), 1, rows2).T(), rows1 * rows2, cols1 * cols2);
    auto tiled_mat2 = af::tile(mat2, rows1, cols1);

    return tiled_mat2 *= resized_mat1;
}

af::array tensor_product_af_opencl(const af::array& mat1, const af::array& mat2)
{
    af::dim4 dims1 = mat1.dims();
    af::dim4 dims2 = mat2.dims();

    af::array out = af::tile(mat2, dims1[0], dims1[1]);
    af::array resized_mat1 = af::resize(dims2[0], dims2[1], mat1, AF_INTERP_LOWER);

    return out *= resized_mat1;
}