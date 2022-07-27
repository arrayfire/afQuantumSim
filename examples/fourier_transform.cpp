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
#include "quantum_gates.h"
#include "quantum_visuals.h"

#include "utils.h"

#include <iostream>

void quantum_fourier_transform_example()
{
    std::cout << "---- Example: Fourier Transform Algorithm ----\n";
    
    uint32_t qubits = 4;
    aqs::QSimulator qs{ qubits };
    aqs::QCircuit qc{ qubits };

    // Represent the value 3 in the fourier transform
    int value = 3;
    for (uint32_t i = 0; i < qubits; ++i)
        qs.qubit(i) = value & (1 << i) ? aqs::QState::one() : aqs::QState::zero();
    qs.generate_statevector();

    std::cout << "Input value: " << value << "\nInput State (Computational Basis):\n";
    aqs::print_statevector(qs);

    aqs::QCircuit fourier_gate = aqs::fourier_transform(qubits);
    qc << aqs::Gate{ fourier_gate, 0 };

    qc.compile();
    qs.simulate(qc);
    std::cout << "Output State (Fourier Basis):\n";
    aqs::print_statevector(qs);

    aqs::QCircuit inverse_fourier_gate = aqs::inverse_fourier_transform(qubits);
    qc.clear();
    qc << aqs::Gate( inverse_fourier_gate, 0 );

    qc.compile();
    qs.simulate(qc);
    std::cout << "Output State after Inverse Fourier (Computational Basis):\n";
    print_statevector(qs);

    std::cout << "Output measurement after Fourier and Inverse Fourier: " << reverse_binary(qs.peek_measure_all(), qubits) << '\n';

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_fourier_transform_example();

    return 0;
}