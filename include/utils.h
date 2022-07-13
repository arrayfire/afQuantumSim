/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#pragma once

#include <arrayfire.h>
#include <string>

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

/**
 * @brief Inserts a given number of bits of a value into another value. It does it per array element
 *
 * @note function only used by gen_index and ControlCircuitGate
 * 
 * @param current to which values the bits will be inserted
 * @param insert what are the values of the bits to be inserted
 * @param position the position where the bits will be inserted
 * @param offset the number of bits that will be inserted
 * @param array_length length of all the arrays passed
 * @return af::array array with the values of the bits inserted
 */
af::array insert_bits(const af::array& current, const af::array& insert, const af::array& position, const af::array& offset, size_t array_length);

/**
 * @brief Generates the indices of the positions in a CircuitGate where the circuit gate entries will be copied to
 * 
 * @note function only used by CircuitGate and ControlCircuitGate
 * 
 * @param gate_qubit_begin starting position of the gate
 * @param gate_qubit_count number of qubits in the gate
 * @param circuit_qubit_count number of qubits in the circuit
 * @param values indices in the gate
 * @param indices index in the circuit
 * @param array_length length of all the arrays passed
 * @return af::array correct indices where the values of the gate should be set
 */
af::array gen_index(uint32_t gate_qubit_begin, uint32_t gate_qubit_count, uint32_t circuit_qubit_count, const af::array& values, const af::array& indices, size_t array_length);
