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

void quantum_entanglement_example()
{
    std::cout << "---- Example: Entanglement ----\n";

    int reps = 1e4;
    aqs::QSimulator qs{4};
    aqs::QCircuit qc{4};

    std::cout << "\n*** Entanglement of 2 qubits ***\n\n";

    // Entagle 2 qubits
    qc << aqs::H{0};
    qc << aqs::CNot{0, 1};

    qc.compile();
    qs.simulate(qc);
    std::cout << "State of entangled q[0] and q[1]:\n";
    aqs::print_statevector(qs);
    //aqs::print_profile(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 3 qubits ***\n\n";

    //Entangle 3 qubits
    qs.generate_statevector();

    qc << aqs::CNot{0, 2};

    qc.compile();
    qs.simulate(qc);
    std::cout << "State of entangled q[0], q[1], and q[2]:\n";
    aqs::print_statevector(qs);
    //aqs::print_statevector(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 4 qubits ***\n\n";

    //Entangle 4 qubits
    qs.generate_statevector();

    qc << aqs::CNot{0, 3};
    qc.compile();
    qs.simulate(qc);

    std::cout << "State of entangled q[0], q[1], q[2], and q[3]:\n";
    aqs::print_statevector(qs);
    //qs.profile_measure_all(reps);

    std::cout << "\n*** Entanglement of 2 pair of qubits ***\n\n";

    //Entangle 2 pair of qubits
    qs.generate_statevector();
    qc.clear();

    qc << aqs::H{0};
    qc << aqs::H{2};
    qc << aqs::CNot{0, 1};
    qc << aqs::CNot{2, 3};

    qc.compile();
    qs.simulate(qc);

    std::cout << "State of entangled (q[0], q[1]), and (q[2], q[3]):\n";
    aqs::print_statevector(qs);
    //qs.profile_measure_all(reps);

    std::cout << "-------------------------\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_entanglement_example();

    return 0;
}