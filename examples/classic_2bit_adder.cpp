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

void quantum_2bit_adder_example()
{
    std::cout << "---- Example: 2-bit adder ----\n";
    std::cout << "Adding number A = q[1]q[0] and B = q[3]q[2] and storing it in A + B = C = q[6]q[5]q[4]:\n\n";

    aqs::QSimulator qs(7);
    aqs::QCircuit qc(7);

    // Add main bits without carry
    qc << aqs::CCNot(0, 2, 5);
    qc << aqs::CCNot(1, 3, 6);

    // Determine carry bits
    qc << aqs::CNot(0, 2);
    qc << aqs::CNot(1, 3);

    // Add carry bits from input
    qc << aqs::CCNot(3, 5, 6);

    // Add carry bits to output
    qc << aqs::CNot(2, 4);
    qc << aqs::CNot(3, 5);

    // Restore original input bits
    qc << aqs::CNot(0, 2);
    qc << aqs::CNot(1, 3);

    auto print_result = [](int val1, int val2, aqs::QSimulator& qs, aqs::QCircuit& qc) {
        qs.qubit(0) = val1 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(1) = val1 & 2 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(2) = val2 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(3) = val2 & 2 ? aqs::QState::one() : aqs::QState::zero();

        qs.generate_statevector();
        std::cout << "Input: " << binary_string(qs.peek_measure_all(), 7) << '\n';

        qc.compile();
        qs.simulate(qc);

        auto result = qs.peek_measure_all();
        std::string r = binary_string(result, 7);

        std::cout << "Output: " << r << '\n';
        std::cout << r[1] << r[0] << " + " << r[3] << r[2] << " = " << r[6] << r[5] << r[4]
                  << "  (" << val1 << " + " << val2 << " = " << reverse_binary(extract_binary(result, 0, 2), 3) << ")\n\n";
    };

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
            print_result(i, j, qs, qc);
    }

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_2bit_adder_example();

    return 0;
}