/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h" // Contains the core simulation functionality from the library
#include "quantum_visuals.h" // Contains functionality for displaying circuits and states

#include <iostream>

int main(int argc, char** argv)
{
    // Initialize library
    aqs::initialize(argc, argv);

    // Create a 2-qubit Quantum Circuit
    aqs::QCircuit qc{ 2 };

    // Add gates to the circuit
    qc << aqs::H{0} << aqs::CX{ 0 , 1 };

    // Compile the circuit
    qc.compile();

    // Create a 2-qubit Simulator with qubits initialize to the |1> state
    aqs::QSimulator qs{ 2 , aqs::QState::one() };

    // Simulate the circuit with the simulator
    qs.simulate(qc);

    // Print the resulting statevector of the simulation
    std::cout << "\nSimulation statevector result:\n";
    aqs::print_statevector(qs);

    // Profile the simulation for 100 simulations
    std::cout << "\nRunning 100 quantum simulation for the same circuit:\n";
    aqs::print_profile(qs.profile_measure_all(100));

    // Print a visual representation of the circuit
    aqs::print_circuit_text_image(qc, qs);

    return 0;
}