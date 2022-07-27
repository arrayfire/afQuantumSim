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

void quantum_teleportation()
{
    std::cout << "---- Example: Quantum Teleportation ----\n";

    aqs::QCircuit qc{3};
    aqs::QSimulator qs(qc.qubit_count());

    float angle = -atan(0.75f);

    // Initial sender qubit state
    qs.qubit(0) = aqs::QState{0.6f, {0.0f, 0.8f}};
    qs.generate_statevector();

    std::cout << "Sender qubit:\n";
    aqs::print_state(qs.qubit(0));
    std::cout << '\n';

    // Create entangled pair
    qc << aqs::H{1} << aqs::CX{1, 2};
    qc << aqs::Barrier{};
    
    // Sender applies some operations
    qc << aqs::CX{0, 1} << aqs::H{0};

    // Sender measures their qubits
    qc.compile();
    qs.simulate(qc);

    auto q1 = qs.measure(0);
    auto q2 = qs.measure(1);

    std::cout << "Sender first qubit measured: " << q1 << "\n";
    std::cout << "Sender second qubit measured: " << q2 << "\n\n";

    // Receiver decodes the qubit received
    aqs::QCircuit qc0{3};
    aqs::QSimulator qs0 = qs;

    if (q2) qc0 << aqs::X{2};

    if (q1) qc0 << aqs::Z{2};

    qc0.compile();
    qs0.simulate(qc0);

    std::cout << "Profiling the qubit received:\n";
    aqs::print_profile(qs.profile_measure(2, 1e4));

    //This is to test if the state is the same as the sender's state
    qc0 << aqs::RotX{2, acos(0.6f) * 2};
    qc0.compile();
    qs.simulate(qc0);

    std::cout << "\nTesting qubit with a RotX gate:\n";
    aqs::print_profile(qs.profile_measure(2, 1e4));
    std::cout << "Expected qubit measurement with gate operation: |0>\n";

    std::cout << "\n-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_teleportation();

    return 0;
}