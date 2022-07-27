/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h"
#include "quantum_algo.h"
#include "quantum_visuals.h"

#include "utils.h"

#include <iostream>
#include <iomanip>

void quantum_constant_addition()
{
    std::cout << "---- Example: Constant Adder ----\n";
    uint32_t qubits = 5;

    int const_val = 8;
    int input_value = 1;

    aqs::QSimulator qs{qubits};
    aqs::QCircuit qc{qubits};

    // Transform from computational to fourier basis
    qc << aqs::Gate(aqs::fourier_transform(qubits), 0, "QFT");

    // Add constant by setting a specific phase rotation
    aqs::QCircuit add_circuit{qubits};
    for (uint32_t i = 0; i < qubits; ++i)
        add_circuit << aqs::Phase(i, const_val * aqs::pi / static_cast<float>(1 << i));
    qc << aqs::Gate{ add_circuit, 0 };

    // Transform back to computational basis
    qc << aqs::Gate{aqs::inverse_fourier_transform(qubits), 0, "QFT†"};
    for (uint32_t i = 0; i < qubits; ++i)
        qs.qubit(i) = (input_value & (1 << i)) ? aqs::QState::one() : aqs::QState::zero();

    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::zero();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();

    std::cout << "\nInput value: " << input_value << " (" << binary_string(input_value, qubits) << ") ; Constant Value: " << const_val << "\n";
    qs.generate_statevector();
    //std::cout << "Input state:\n";
    //aqs::print_statevector(qs);

    qc.compile();
    qs.simulate(qc);

    //std::cout << "Output state:\n";
    //aqs::print_statevector(qs);
    //aqs::print_profile(qs.profile_measure_all(1e4));

    std::cout << "Output value: " << binary_string(reverse_binary(qs.peek_measure_all(), qubits), qubits) << "\n";
    std::cout << "Result: " << input_value << " + " << const_val << " = " << reverse_binary(qs.peek_measure_all(), qubits) << "\n\n";

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

void quantum_two_addition()
{
    std::cout << "---- Example: Two Number Adder ----\n";

    uint32_t input = 4;

    aqs::QCircuit qc{2 * input + 1};
    aqs::QSimulator qs{2 * input + 1};

    // Transform to fourier basis
    // Execute Fourier Transform on b input
    qc << aqs::Gate{ aqs::fourier_transform(input + 1), input, "QFT" };

    // Execute Fourier Transform on a input stored on b input
    // This is equivalent to add phase rotations
    for (uint32_t i = 0; i < input; ++i)
        qc << aqs::CPhase{ i, input * 2, aqs::pi / (1 << (input - i)) };

    for (int32_t i = input; i >= 0; --i)
    {
        for (uint32_t j = 0; j < i; ++j)
            qc << aqs::CPhase{ j, input + i - 1, aqs::pi / (1 << (i - j - 1)) };
    }

    // Return back to computational basis
    qc << aqs::Gate{ aqs::inverse_fourier_transform(input + 1), input, "QFT†" };

    std::cout << "\nAddition of two numbers of " << input << " bits\n"
              << "First " << input << " qubits are number A and the next " << input << " qubits are number B\n"
              << "Algorithm computes C = A + B, and stores C in the last " << input + 1 << " qubits\n\n";

    for (int i = 0; i < fast_pow2(input * 2); ++i)
    {
        int valA = i % fast_pow2(input);
        int valB = i / fast_pow2(input);

        // Setup the registers with the values to add
        for (uint32_t j = 0; j < input; ++j)
            qs.qubit(j) = (valA & (1 << j)) ? aqs::QState::one() : aqs::QState::zero();
        for (uint32_t j = 0; j < input; ++j)
            qs.qubit(j + input) = (valB & (1 << j)) ? aqs::QState::one() : aqs::QState::zero();

        // Reset the carry bit
        qs.qubit(2 * input) = aqs::QState::zero();

        // Prepare the simulation statevector
        qs.generate_statevector();

        // Compile and simulate
        qc.compile();
        qs.simulate(qc);

        // Select the last input + 1 bits and flip them to obtain the result in binary
        int result = reverse_binary(extract_binary(qs.peek_measure_all(), 0, input), input + 1);

        std::string out = binary_string(qs.peek_measure_all(), input + 1);
        std::reverse(out.begin(), out.end());
        
        std::cout << binary_string(valA, input) << " + " << binary_string(valB, input) << " = " << out
                  << " ( " << std::setw(2) << valA << " + " << std::setw(2) << valB << " = " << std::setw(2) << result << " )\n";
    }

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "\n-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_constant_addition();
    quantum_two_addition();   

    return 0;
}