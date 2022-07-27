/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h"
#include "quantum_visuals.h"

#include <iostream>

void quantum_classic_xor_example()
{
    std::cout << "---- Example: XOR Gate ----\n";
    std::cout << "XOR qubits 1 and 2 with qubit 0:\n\n";

    int reps = 1e4;

    aqs::QSimulator qs(3);
    aqs::QCircuit qc(3);
    
    qc << aqs::Xor(0, 1);
    qc << aqs::Xor(0, 2);

    // XOR with state 0
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_statevector();

    std::cout << "XOR when q[0] = 0\nInitial State:\n";
    print_statevector(qs);

    qc.compile();
    qs.simulate(qc);

    std::cout << "Output State:\n";
    print_statevector(qs);
    std::cout << "Expected: 001 -> 001\n";
    std::cout << "Measured: 001 -> " << binary_string(qs.peek_measure_all(), 3) << "\n"; 
    //qs.profile_measure_all(reps);

    // XOR with state 1
    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_statevector();

    std::cout << "\n*******\n\n";
    std::cout << "XOR when q[0] = 1\nInitial State:\n";
    print_statevector(qs);

    qs.simulate(qc);

    std::cout << "Output State:\n";
    print_statevector(qs);
    std::cout << "Expected: 101 -> 110\n";
    std::cout << "Measured: 101 -> " << binary_string(qs.peek_measure_all(), 3) << "\n";
    //qs.profile_measure_all(reps);

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

void quantum_classic_or_example()
{
    std::cout << "---- Example: OR Gate ----\n";
    std::cout << "OR qubits 0 and 1 and output it into qubit 3,\nOR qubits 0 and 2 and output it into qubit 4:\n\n";
    int reps = 1e4;

    aqs::QSimulator qs(5);
    aqs::QCircuit qc(5);
    
    qc << aqs::Or(0, 1, 3);
    qc << aqs::Or(0, 2, 4);

    // OR with state 0
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_statevector();
    std::cout << "OR when q[0] = 0 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_statevector(qc);

    qc.compile();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_statevector(qs);
    std::cout << "Expected: 00100 -> 00101\n";
    std::cout << "Measured: 00100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n"; 
    //qc.profile_measure_all(reps);

    std::cout << "\n******\n\n";

    // OR with state 1
    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_statevector();
    std::cout << "OR when q[0] = 1 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_statevector(qc);

    qc.compile();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_statevector(qs);
    std::cout << "Expected: 10100 -> 10111\n";
    std::cout << "Measured: 10100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n";
    //qs.profile_measure_all(reps);

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

void quantum_classic_and_example()
{
    std::cout << "---- Example: AND Gate ----\n";
    std::cout << "AND qubits 0 and 1 and output it into qubit 3,\nOR qubits 0 and 2 and output it into qubit 4:\n\n";

    aqs::QCircuit qc{ 5 };
    aqs::QSimulator qs{ 5 };

    qc << aqs::And(0, 1, 3);
    qc << aqs::And(0, 2, 4);

    // OR with state 0
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_statevector();
    std::cout << "AND when q[0] = 0 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_statevector(qc);

    qc.compile();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_statevector(qs);
    std::cout << "Expected: 00100 -> 00100\n";
    std::cout << "Measured: 00100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n"; 
    //qc.profile_measure_all(reps);

    std::cout << "\n******\n\n";

    // OR with state 1
    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_statevector();
    std::cout << "AND when q[0] = 1 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_statevector(qc);

    qc.compile();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_statevector(qs);
    std::cout << "Expected: 10100 -> 10101\n";
    std::cout << "Measured: 10100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n";
    //qs.profile_measure_all(reps);

    //aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_classic_or_example();

    quantum_classic_xor_example();

    quantum_classic_and_example();

    return 0;
}