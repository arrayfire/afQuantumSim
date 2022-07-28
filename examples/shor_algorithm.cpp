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
#include <iomanip>
#include <iostream>
#include <set>

void quantum_shor_algorithm()
{
    std::cout << "---- Example: Shor's Algorithm ----\n";

    uint32_t output = 8;

    aqs::QCircuit qc{ output + 4 };
    aqs::QSimulator qs{ output + 4 };

    qs.qubit(qs.qubit_count() - 1) = aqs::QState::one();

    std::cout << "Adding Hadamard Gates...\n";

    aqs::QCircuit h1(output / 2), h2(output - (output / 2));
    for (uint32_t i = 0; i < h1.qubit_count(); ++i)
        h1 << aqs::H{i};
    for (uint32_t i = 0; i < h2.qubit_count(); ++i)
        h2 << aqs::H{i};
    qc << aqs::Gate{ h1, 0 };
    qc << aqs::Gate{ h2, h1.qubit_count() };

    // Create a^(2^j) mod 15 gate
    auto a2jmod15 = [](int a, int j) {
        aqs::QCircuit temp(4);
        if (a == 13 || a == 2)
        {
            temp << aqs::Swap(0, 1);
            temp << aqs::Swap(1, 2);
            temp << aqs::Swap(2, 3);
        }
        if (a == 7 || a == 8)
        {
            temp << aqs::Swap(2, 3);
            temp << aqs::Swap(1, 2);
            temp << aqs::Swap(0, 1);
        }
        if (a == 11 || a == 4)
        {
            temp << aqs::Swap(1, 3);
            temp << aqs::Swap(0, 2);
        }
        if (a == 7 || a == 11 || a == 13)
        {
            for (uint32_t j = 0; j < 4; ++j)
                temp << aqs::X{j};
        }

        for (uint32_t i = 0; i < j; ++i)
            temp << aqs::Gate{ temp, 0 };

        return temp;
    };

    // Biggest number than can be represented with 4 bits
    int N = 15;

    // Factor integer `a`
    int a = 13;

    for (uint32_t i = 0; i < output; ++i)
    {
        std::cout << "Adding Control_mod15 a^(2^" << i << ") Gate...\n";
        qc << aqs::ControlGate{ a2jmod15(a, i), output - 1 - i, output, "a^(2^" + std::to_string(i) + ") mod 15" };
    }

    std::cout << "Adding Fourier Transform Gate...\n";
    qc << aqs::Gate(aqs::inverse_fourier_transform(output), 0, "QFT†");

    std::cout << "Generating Statevector...\n";
    qs.generate_statevector();

    std::cout << "Generating circuit...\n";
    qc.compile();

    std::cout << "Executing Simulation...\n";
    qs.simulate(qc);

    std::cout << "Profiling...\n";
    int reps = 1e4;
    auto profile = qs.profile_measure_all(reps);

    std::vector<float> qubit_probs(output);
    for (int i = 0; i < qubit_probs.size(); ++i)
        qubit_probs[i] = qs.qubit_probability_true(i);

    // Measure all the output qubits
    std::vector<std::pair<int, float>> state_probs(fast_pow2(output), {0 , 1.0f});
    for (int i = 0; i < state_probs.size(); ++i)
    {
        state_probs[i].first = i;
        for (int j = 0; j < output; ++j)
            state_probs[i].second *= (i & (1 << j)) ? qubit_probs[j] : (1.0f - qubit_probs[j]);
    }

    // Sort the states from the measured qubits in descending probability
    std::sort(state_probs.begin(), state_probs.end(), [](const std::pair<int, float>& a, const std::pair<int, float>& b){
        return b.second < a.second;
    });

    std::set<int> r_values;

    // Find the most probable values for r
    int top_count = 4;
    std::cout << "\nTop " << top_count << " Results:\n";
    for (int i = 0; i < top_count; ++i)
    {
        auto state = state_probs[i].first;
        auto prob = state_probs[i].second;
        auto base = fast_pow2(output);
        auto phase = static_cast<double>(state_probs[i].first) / static_cast<double>(base);
        auto frac = approximate_fraction(phase, 15);

        // Store all unique values for r
        r_values.insert(static_cast<int>(frac.second));

        std::cout << binary_string(state, output) << " = " << std::setw(4) << state << " (Prob: " << std::setprecision(3) << prob
                  << std::setprecision(7) << ") - Phase: " << std::setw(4) << std::right << state << " / " << base << " = "
                  << std::setw(4) << phase << " = " << frac.first << " / " << frac.second << "\n";
    }

    // Execute factoring with the values of r
    int guesses = 1;
    std::cout << "\nValue of a = " << a << "\n\n";
    for (const auto& r : r_values)
    {
        std::cout << "Attempt #" << guesses++ << "\n r = " << r << "\n"; 
        if (r % 2)
        {
            std::cout << "Found r is odd, continuing next attempt...\n\n";
            continue;
        }

        int guess1 = static_cast<int>(gcd(pow(a, r / 2) - 1, N));
        int guess2 = static_cast<int>(gcd(pow(a, r / 2) + 1, N));
        std::cout << "Guessed factors: " << guess1 << " and " << guess2 << "\n";

        if (guess1 != 1 && guess1 != N)
            std::cout << "**Non-trivial factor: " << guess1 << "**\n";
        if (guess2 != 1 && guess2 != N)
            std::cout << "**Non-trivial factor: " << guess2 << "**\n";
        std::cout << "\n";
    }

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    quantum_shor_algorithm();

    return 0;
}