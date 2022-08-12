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

#include <algorithm>
#include <iostream>
#include <iomanip>

void quantum_grover_example()
{
    std::cout << "---- Example: Grover Search Algorithm ----\n";

    uint32_t qubits = 4;
    int search_space = fast_pow2(qubits);
    int solution_count = 1;
    int marked_solution = 5;
    int iterations = aqs::pi * std::sqrt(static_cast<float>(search_space) / static_cast<float>(solution_count)) / 4.f;
    
    aqs::QSimulator qs{ qubits };
    aqs::QCircuit qc{ qubits };

    // Create an oracle that marks a solution
    aqs::QCircuit oracle_gate = aqs::grover_oracle(qubits, marked_solution);

    // Create a grover search gate with the create oracle
    aqs::QCircuit grover_search_gate = aqs::grover_search(qubits, oracle_gate, iterations, "Oracle");

    qc << aqs::Gate{ grover_search_gate, 0 };

    qc.compile();
    qs.simulate(qc);

    std::cout << "Searching for <" << marked_solution << "> in " << search_space << " numbers using " << iterations
              << " iterations of Grover's Algorithm\n\n";
    //std::cout << "Grover Output State:\n";
    //print_statevector(qs);

    int reps = 1e4;
    std::cout << "Profile of " << reps << " measurements:\n";
    auto profile = qs.profile_measure_all(reps);
    //aqs::print_profile(profile);

    std::vector<std::pair<int, int>> results(profile.size());
    for (int i = 0; i < profile.size(); ++i)
        results[i] = { i , profile[i] };

    // Find the most probable state
    std::sort(results.begin(), results.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return b.second < a.second;
    });

    // Output the state that where amplified
    int top_result_count = 4;
    std::cout << "Top " << top_result_count << " Results:\n";
    for (int i = 0; i < top_result_count; ++i)
        std::cout << binary_string(results[i].first, qubits) << " (" << std::setw(1 + qubits * 3 / 10) <<  reverse_binary(results[i].first, qubits) << ")"
                  << " : " << std::setw(4) << static_cast<float>(results[i].second * 100) / static_cast<float>(reps) << "%\n";

    std::cout << "\n**Search found <" << reverse_binary(results[0].first, qubits) << "> as the solution**\n";

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_grover_example();

    return 0;
}