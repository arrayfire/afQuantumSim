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

#include <algorithm>
#include <iostream>

void quantum_counting()
{
    std::cout << "---- Example: Quantum Grover Counting ----\n";

    uint32_t output_count = 4;
    uint32_t search_count = 4;

    aqs::QCircuit qc{ output_count + search_count };
    aqs::QSimulator qs{ output_count + search_count };

    for (uint32_t i = 0; i < output_count + search_count; ++i)
        qc << aqs::H{i};

    aqs::QCircuit grover_iter{ 4 };

    // Add Oracle
    grover_iter << aqs::H{2} << aqs::H{3} << aqs::Gate{aqs::NControl_Gate(4, {0,1}, 2, aqs::X::gate()), 0} << aqs::H{2} << aqs::X{2}
                << aqs::Gate{aqs::NControl_Gate(4, {0, 2}, 3, aqs::X::gate()), 0 } << aqs::X{2} << aqs::H{3} << aqs::X{1}
                << aqs::X{3} << aqs::H{2}
                << aqs::Gate{aqs::NControl_Gate(4, {0,1,3}, 2, aqs::X::gate()), 0} << aqs::X{1} << aqs::X{3} << aqs::H{2};

    // Add Amplifier
    grover_iter << std::vector<aqs::H>{0,1,2} << std::vector<aqs::X>{0,1,2} << aqs::Z{3}
                << aqs::Gate{ aqs::NControl_Gate(4, {0,1,2}, 3, aqs::X::gate()), 0}
                << std::vector<aqs::X>{0,1,2} << std::vector<aqs::H>{0,1,2} << aqs::Z{3};

    // Add Control Grover Iterations
    for (int32_t i = output_count - 1; i >= 0; --i)
    {
        qc << aqs::ControlGate{ grover_iter, output_count - 1 - i, output_count, "G" };
        grover_iter << aqs::Gate{ grover_iter, 0 };
    }

    // Execute phase estimation
    qc << aqs::Gate{ aqs::inverse_fourier_transform(output_count), 0, "QFT†" };
    qc.compile();

    qs.simulate(qc);
    print_circuit_text_image(qc, qs);

    std::vector<float> qubit_probs(output_count);
    for (int i = 0; i < qubit_probs.size(); ++i)
        qubit_probs[i] = qs.qubit_probability_true(i);

    // Find the probabilities for the states measured from output count
    std::vector<std::pair<uint32_t, float>> state_probs(fast_pow2(output_count), {0 , 1.0f});
    for (int i = 0; i < state_probs.size(); ++i)
    {
        state_probs[i].first = i;
        for (int j = 0; j < output_count; ++j)
            state_probs[i].second *= (i & (1 << j)) ? qubit_probs[j] : (1.0f - qubit_probs[j]);
    }

    // Sort the state probabilities descending
    std::sort(state_probs.begin(), state_probs.end(), [](const std::pair<uint32_t, float>& a, const std::pair<uint32_t, float>& b){
        return b.second < a.second;
    });

    // Find the most likely state
    auto likely_value = state_probs[0].first;
    auto phase = (float)likely_value * aqs::pi * 2.f / (float)(1 << output_count);

    auto search_space = 1 << search_count;
    auto solution_count = (1 << search_count) * std::cos(phase / 2.f) * std::cos(phase / 2.f);
    std::cout << "Likely output from phase estimation: " << likely_value << "(" << binary_string(likely_value, output_count) << ")\n";
    std::cout << "Estimated phase: " << phase << '\n';
    std::cout << "Estimated number of solutions: " << solution_count << '\n';

    auto error = (std::sqrt(2.f * search_space * (search_space - solution_count)) +
                (float)search_space / (float)(1 << output_count)) / (float)(1 << (output_count - 1));
    std::cout << "Error: " << error << "\n\n";

    std::cout << "Expected number of solutions: 5\n";

    std::cout << "\n-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_counting();

    return 0;
}