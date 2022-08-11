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

#include <iostream>
#include <chrono>
#include <set>
#include <algorithm>

void qft(uint32_t qubits)
{
    aqs::QCircuit qc = aqs::fourier_transform(qubits);
    aqs::QSimulator qs{ qubits };

    //qc.compile();
    qs.simulate(qc);

    auto result = qs.statevector();
    result.eval();
}

void shor(uint32_t qubits)
{
    uint32_t output = 8;

    aqs::QCircuit qc{ output + 4 };
    aqs::QSimulator qs{ output + 4 };

    qs.qubit(qs.qubit_count() - 1) = aqs::QState::one();

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
        qc << aqs::ControlGate{ a2jmod15(a, i), output - 1 - i, output, "a^(2^" + std::to_string(i) + ") mod 15" };
    }

    qc << aqs::Gate(aqs::inverse_fourier_transform(output), 0, "QFT†");

    //qc.compile();
    qs.simulate(qc);

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

    // Find the top probable states
    std::size_t top_count = 4;
    for (std::size_t i = 0; i < top_count; ++i)
    {
        std::size_t max_index = i;
        float max_value = state_probs[i].second;
        for (std::size_t j = i + 1; j < state_probs.size(); ++j)
        {
            if (state_probs[j].second > max_value)
            {
                max_value = state_probs[j].second;
                max_index = j;
            }
        }
        std::swap(state_probs[max_index], state_probs[i]);
    }

    std::set<int> r_values;

    // Find the most probable values for r
    for (int i = 0; i < top_count; ++i)
    {
        auto state = state_probs[i].first;
        auto prob = state_probs[i].second;
        auto base = fast_pow2(output);
        auto phase = static_cast<double>(state_probs[i].first) / static_cast<double>(base);
        auto frac = approximate_fraction(phase, 15);

        // Store all unique values for r
        r_values.insert(static_cast<int>(frac.second));
    }

    // Execute factoring with the values of r
    int guesses = 1;
    
    for (const auto& r : r_values)
    {
        if (r % 2)
            continue;

        int guess1 = static_cast<int>(gcd(pow(a, r / 2) - 1, N));
        int guess2 = static_cast<int>(gcd(pow(a, r / 2) + 1, N));
    }
}

void grover(uint32_t qubits)
{
    uint32_t marked_state = 5;
    uint32_t solution_count = 1;
    uint32_t search_space = 1 << qubits;

    uint32_t iterations = static_cast<uint32_t>(aqs::pi * std::sqrt(search_space / solution_count) / 4.f);
    aqs::QCircuit qc = aqs::grover_search(qubits, aqs::grover_oracle(qubits, marked_state), iterations);
    aqs::QSimulator qs{ qubits };

    //qc.compile();
    qs.simulate(qc);

    auto result = qs.statevector();
    result.eval();
}

void benchmark(uint32_t qubits, uint32_t count, void (*func)(uint32_t), const char* name)
{
    int64_t total = 0;
    int64_t total_squared = 0;

    for (uint32_t i = 0; i < count; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();
        func(qubits);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        total += duration;
        total_squared += duration * duration;
    }

    float average = total / (float)count;
    float squared_average = total_squared / (float)count;

    std::cout << name << "\n";
    std::cout << "Total time: " << total / 1000.f << " ms; Average Time: " << average / 1000.f
              << " ms; Std. Dev: " << count / (float)(count - 1) * std::sqrt(squared_average - average * average) / 1000.f << " ms; Count: " << count << '\n' << std::endl;
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);

    const uint32_t qubits = 10;
    const uint32_t count = 100;

    aqs::clear_circuit_cache();
    benchmark(qubits, count, grover, "Grover");

    aqs::clear_circuit_cache();
    benchmark(qubits, count, qft, "QFT");

    aqs::clear_circuit_cache();
    benchmark(qubits, count, shor, "Shor");

    return 0;
}