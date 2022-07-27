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

#include "utils.h"

#include <iostream>

void quantum_superposition_example()
{
    std::cout << "---- Example: Superposition ----\n";

    int reps = 1e4;
    aqs::QSimulator qs{4};
    aqs::QCircuit qc{4};

    // All qubits in definite |0> state
    qs.generate_statevector();
    std::cout << "Statevector with all qubits in |0> state:\n";
    print_statevector(qs);
    //aqs::print_profile(qs.profile_measure_all(reps));

    // Put qubits 0 and 2 into superposition
    qc << aqs::H{0};
    qc << aqs::H{2};

    qc.compile();
    qs.simulate(qc);

    std::cout << "Statevector with q[0] and q[2] in superposition using Hadamard:\n";
    print_statevector(qs);
    //print_profile_measure(qs.profile_measure_all(reps));

    // Put qubit 1 into a superposed state and qubit 3 into definite |1> state
    qs.qubit(1) = aqs::QState{{0.0f, 0.6f} , { 0.8f, 0.0f }};
    qs.qubit(3) = aqs::QState::one();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    std::cout << "Statevector with q[0] and q[2] in superposition using Hadamard, q[1] = 0.6i|0> + 0.8|1>, and q[3] = |1>:\n";
    print_statevector(qs);
    //print_profile(qs.profile_measure_all(reps));

    std::cout << "-------------------------\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_superposition_example();

    return 0;
}