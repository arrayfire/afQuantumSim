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

/**
 * @brief Returns a string from the repetitions of the passed string
 *
 * @param n number of repetitions
 * @param str string to repeat
 * @return std::string
 */
std::string repeat(std::size_t n, const std::string& str);

/**
 * @brief Returns the length of a UTF8 string stored in the string
 *
 * @param str string to obtain the length of
 * @return std::size_t
 */
std::size_t utf8str_len(std::string str);

/**
 * @brief Executes the tensor product between two arrayfire arrays
 *
 * @param rhs lhs array
 * @param rhs rhs array
 * @return af::array result dims[0] = lhs.dims[0] * rhs.dims[1], dims[1] =
 * lhs.dims[1] * rhs.dims[1]
 */
af::array tensor_product(const af::array& lhs, const af::array& rhs);

/**
 * @brief Returns the string representation of the given measurement
 *
 * @param val value of the measurement
 * @param length the number of qubits/the length of the bits to read (should be
 * less than or equal to 32)
 * @return std::string
 */
std::string binary_string(uint32_t val, int length);

/**
 * @brief Reverses the bits in the given integer value
 *
 * @param val value to reverse
 * @param length number of bits to be reverse (less than or equal to 32)
 * @return uint32_t result
 */
uint32_t reverse_binary(uint32_t val, int length) noexcept;

/**
 * @brief Returns the number formed from the range of bits selected
 *
 * @param val value from which the number will be extracted
 * @param first index of the first bit
 * @param last index of the last bit (should be greater than first)
 * @return uint32_t result
 */
uint32_t extract_binary(uint32_t val, int first, int last) noexcept;

/**
 * @brief Al
 *
 * @param value
 * @param max_denominator
 * @return std::pair<int64_t, int64_t>
 */
std::pair<int64_t, int64_t> approximate_fraction(double value,
                                                 int64_t max_denominator);

/**
 * @brief Returns the Greates Common Divisor (GCD) of two numbers
 *
 * @details Uses euclid's algorithm
 *
 * @param a
 * @param b
 * @return int64_t
 */
int64_t gcd(int64_t a, int64_t b);

/**
 * @brief Returns the power of 2 that correspond to the given exponent
 *
 * @param pow Power/Exponent
 * @return uint32_t
 */
static inline uint32_t fast_pow2(uint32_t pow) { return 1 << pow; }

/**
 * @brief Returns the truncated integer logarithm base 2 of the given number
 *
 * @param val number to log
 * @return uint32_t
 */
static inline uint32_t fast_log2(uint32_t val) {
    int counter = -(val != 0);
    for (; val; counter++) val >>= 1;
    return counter;
}

/**
 * @brief Inserts a given number of bits of a value into another value. It does
 * it per array element
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
af::array insert_bits(const af::array& current, const af::array& insert,
                      const af::array& position, const af::array& offset,
                      size_t array_length);

/**
 * @brief Generates the indices of the positions in a CircuitGate where the
 * circuit gate entries will be copied to
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
af::array gen_index(uint32_t gate_qubit_begin, uint32_t gate_qubit_count,
                    uint32_t circuit_qubit_count, const af::array& values,
                    const af::array& indices, size_t array_length);
