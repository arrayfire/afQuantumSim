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

#include <algorithm>
#include <iostream>

void quantum_phase_estimation_example()
{
    std::cout << "---- Example: Phase Estimation Algorithm ----\n";

    uint32_t qcount = 7;
    aqs::QSimulator qs(qcount + 1);
    aqs::QCircuit qc(qcount + 1);

    qs.qubit(qcount) = aqs::QState::one();
    qs.generate_statevector();

    for (uint32_t i = 0; i < qcount; ++i)
        qc << aqs::H(i);

    float theta = 1.f/5.f;
    float angle = aqs::pi * 2.f * theta;
    
    std::cout << "\nEstimating a phase shift of " << theta << '\n';

    for (uint32_t i = 0; i < qcount; ++i)
        qc << aqs::CPhase(i, qcount, angle * static_cast<float>(1 << (qcount - 1 - i)));

    qc.compile();
    qs.simulate(qc);
    //std::cout << "Phase shifted state:\n";
    //aqs::print_statevector(qs);
    qs.generate_statevector();

    qc << aqs::Gate(aqs::inverse_fourier_transform(qcount), 0, "QFT†");

    qc.compile();
    qs.simulate(qc);
    //std::cout << "Phase estimation state output: ";
    //aqs::print_statevector(qs);

    auto measurements = qs.profile_measure_all(1e4);
    auto common = std::max_element(measurements.begin(), measurements.end());
    int likely_state = std::distance(measurements.begin(), common);
    std::cout << "Most likely state: " << binary_string(likely_state, qcount + 1)
              << "\t(The reverse of the first " << qcount << " bits is the estimated numerator)\n";

    int estimated_numerator = reverse_binary(likely_state & ~1, qcount + 1);
    std::cout << "Estimated numerator: " << estimated_numerator << "\n";
    std::cout << "Estimated phase shift: " << static_cast<float>(estimated_numerator) / static_cast<float>(1 << qcount)
              << " (" << estimated_numerator << " / " << (1 << qcount) << ")\n\n";

    
    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;
 
    quantum_phase_estimation_example();

    return 0;
}