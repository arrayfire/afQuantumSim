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

#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>

void quantum_classic_xor_example()
{
    using namespace aqs;
    std::cout << "---- Example: XOR Gate ----\n";
    std::cout << "XOR qubits 1 and 2 with qubit 0:\n\n";

    int reps = 1e4;

    aqs::QSimulator qs(3);
    aqs::QCircuit qc(3);
    
    qc << Xor(0, 1);
    qc << Xor(0, 2);

    // XOR with state 0
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_global_state();

    std::cout << "XOR when q[0] = 0\nInitial State:\n";
    print_global_state(qs);

    qc.generate_circuit();
    qs.simulate(qc);

    std::cout << "Output State:\n";
    print_global_state(qs);
    std::cout << "Expected: 001 -> 001\n";
    std::cout << "Measured: 001 -> " << binary_string(qs.peek_measure_all(), 3) << "\n"; 
    //qs.profile_measure_all(reps);

    // XOR with state 1
    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_global_state();

    std::cout << "\n*******\n\n";
    std::cout << "XOR when q[0] = 1\nInitial State:\n";
    print_global_state(qs);

    qs.simulate(qc);

    std::cout << "Output State:\n";
    print_global_state(qs);
    std::cout << "Expected: 101 -> 110\n";
    std::cout << "Measured: 101 -> " << binary_string(qs.peek_measure_all(), 3) << "\n";
    //qs.profile_measure_all(reps);

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

void quantum_classic_or_example()
{
    using namespace aqs;
    std::cout << "---- Example: OR Gate ----\n";
    std::cout << "OR qubits 0 and 1 and output it into qubit 3,\nOR qubits 0 and 2 and output it into qubit 4:\n\n";
    int reps = 1e4;

    aqs::QSimulator qs(5);
    aqs::QCircuit qc(5);
    
    qc << Or(0, 1, 3);
    qc << Or(0, 2, 4);

    // OR with state 0
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_global_state();
    std::cout << "OR when q[0] = 0 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_global_state(qc);

    qc.generate_circuit();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_global_state(qs);
    std::cout << "Expected: 00100 -> 00101\n";
    std::cout << "Measured: 00100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n"; 
    //qc.profile_measure_all(reps);

    std::cout << "\n******\n\n";

    // OR with state 1
    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();
    qs.generate_global_state();
    std::cout << "OR when q[0] = 1 , q[1] = 0 , q[2] = 1\n";
    //std::cout << "Initial State:\n";
    //print_global_state(qc);

    qc.generate_circuit();
    qs.simulate(qc);

    //std::cout << "Output State:\n";
    //print_global_state(qs);
    std::cout << "Expected: 10100 -> 10111\n";
    std::cout << "Measured: 10100 -> " << binary_string(qs.peek_measure_all(), 5) << "\n";
    //qs.profile_measure_all(reps);

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

void quantum_2bit_adder_example()
{
    using namespace aqs;
    std::cout << "---- Example: 2-bit adder ----\n";
    std::cout << "Adding number A = q[1]q[0] and B = q[3]q[2] and storing it in A + B = C = q[6]q[5]q[4]:\n\n";
    aqs::QSimulator qs(7);
    aqs::QCircuit qc(7);

    // Add main bits without carry
    qc << CCNot(0, 2, 5);
    qc << CCNot(1, 3, 6);

    // Determine carry bits
    qc << CNot(0, 2);
    qc << CNot(1, 3);

    // Add carry bits from input
    qc << CCNot(3, 5, 6);

    // Add carry bits to output
    qc << CNot(2, 4);
    qc << CNot(3, 5);

    // Restore original input bits
    qc << CNot(0, 2);
    qc << CNot(1, 3);

    auto print_result = [](int val1, int val2, aqs::QSimulator& qs, aqs::QCircuit& qc) {
        qs.qubit(0) = val1 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(1) = val1 & 2 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(2) = val2 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(3) = val2 & 2 ? aqs::QState::one() : aqs::QState::zero();

        qs.generate_global_state();
        std::cout << "Input: " << binary_string(qs.peek_measure_all(), 7) << '\n';

        qc.generate_circuit();
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

    std::cout << qc.representation() << "\n";
    //aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n";
}

void quantum_entanglement_example()
{
    using namespace aqs;
    std::cout << "---- Example: Entanglement ----\n";

    int reps = 1e4;
    aqs::QSimulator qs(4);
    aqs::QCircuit qc(4);

    std::cout << "\n*** Entanglement of 2 qubits ***\n\n";
    // Entagle 2 qubits
    qc << H(0);
    qc << CNot(0, 1);
    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "State of entangled q[0] and q[1]:\n";
    aqs::print_global_state(qs);
    //aqs::print_profile(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 3 qubits ***\n\n";
    //Entangle 3 qubits
    qc.generate_circuit();
    qs.generate_global_state();
    qc << CNot(0, 2);
    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "State of entangled q[0], q[1], and q[2]:\n";
    aqs::print_global_state(qs);
    //aqs::print_global_state(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 4 qubits ***\n\n";
    //Entangle 4 qubits
    qs.generate_global_state();
    qc << CNot(0, 3);
    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "State of entangled q[0], q[1], q[2], and q[3]:\n";
    aqs::print_global_state(qs);
    //qs.profile_measure_all(reps);

    std::cout << "\n*** Entanglement of 2 pair of qubits ***\n\n";
    //Entangle 2 pair of qubits
    qs.generate_global_state();
    qc.reset_circuit();
    qc << H(0);
    qc << H(2);
    qc << CNot(0, 1);
    qc << CNot(2, 3);
    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "State of entangled (q[0], q[1]), and (q[2], q[3]):\n";
    aqs::print_global_state(qs);
    //qs.profile_measure_all(reps);

    std::cout << "-------------------------\n";
}

void quantum_superposition_example()
{
    using namespace aqs;
    std::cout << "---- Example: Superposition ----\n";

    int reps = 1e4;
    aqs::QSimulator qs(4);
    aqs::QCircuit qc(4);

    // All qubits in definite |0> state
    qs.generate_global_state();
    std::cout << "Global State with all qubits in |0> state:\n";
    print_global_state(qs);
    //aqs::print_profile(qs.profile_measure_all(reps));

    // Put qubits 0 and 2 into superposition
    qc << H(0);
    qc << H(2);
    qc.generate_circuit();
    qs.simulate(qc);

    std::cout << "Global State with q[0] and q[2] in superposition using Hadamard:\n";
    print_global_state(qs);
    //print_profile_measure(qs.profile_measure_all(reps));

    // Put qubit 1 into a superposed state and qubit 3 into definite |1> state
    qs.qubit(1) = aqs::QState{{0.0f, 0.6f} , { 0.8f, 0.0f }};
    qs.qubit(3) = aqs::QState::one();
    qs.generate_global_state();
    qc.generate_circuit();
    qs.simulate(qc);

    std::cout << "Global State with q[0] and q[2] in superposition using Hadamard, q[1] = 0.6i|0> + 0.8|1>, and q[3] = |1>:\n";
    print_global_state(qs);
    //print_profile(qs.profile_measure_all(reps));

    std::cout << "-------------------------\n";
}

void quantum_grover_example()
{
    using namespace aqs;
    std::cout << "---- Example: Grover Search Algorithm ----\n";

    int qubits = 4;
    int search_space = fast_pow2(qubits);
    int solution_count = 1;
    int marked_solution = 5;
    int iterations = aqs::pi * std::sqrt(static_cast<float>(search_space) / static_cast<float>(solution_count)) / 4.f;
    
    aqs::QSimulator qs(qubits);
    aqs::QCircuit qc(qubits);

    aqs::QCircuit oracle_gate = aqs::grover_oracle(qubits, marked_solution);
    aqs::QCircuit grover_search_gate = aqs::grover_search(qubits, oracle_gate, iterations);
    std::cout << grover_search_gate.representation() << std::endl;
    qc << Gate(grover_search_gate, 0);

    qc.generate_circuit();
    qs.simulate(qc);

    std::cout << "Searching for <" << marked_solution << "> in " << search_space << " numbers using " << iterations
              << " iterations of Grover's Algorithm\n\n";
    //std::cout << "Grover Output State:\n";
    //print_global_state(qs);

    int reps = 1e4;
    std::cout << "Profile of " << reps << " measurements:\n";
    auto profile = qs.profile_measure_all(reps);
    //aqs::print_profile(profile);

    std::vector<std::pair<int, int>> results(profile.size());
    for (int i = 0; i < profile.size(); ++i)
        results[i] = { i , profile[i] };

    std::sort(results.begin(), results.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return b.second < a.second;
    });

    int top_result_count = 4;
    std::cout << "Top " << top_result_count << " Results:\n";
    for (int i = 0; i < top_result_count; ++i)
        std::cout << binary_string(results[i].first, qubits) << " (" << std::setw(1 + qubits * 3 / 10) <<  reverse_binary(results[i].first, qubits) << ")"
                  << " : " << std::setw(4) << static_cast<float>(results[i].second * 100) / static_cast<float>(reps) << "%\n";

    std::cout << "\n**Search found <" << reverse_binary(results[0].first, qubits) << "> as the solution**\n";

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

void quantum_fourier_transform_example()
{
    using namespace aqs;
    std::cout << "---- Example: Fourier Transform Algorithm ----\n";
    int qubits = 4;
    aqs::QSimulator qs(qubits);
    aqs::QCircuit qc(qubits);

    int value = 3;
    for (int i = 0; i < qubits; ++i)
        qs.qubit(i) = value & (1 << i) ? aqs::QState::one() : aqs::QState::zero();
    qs.generate_global_state();

    std::cout << "Input value: " << value << "\nInput State (Computational Basis):\n";
    aqs::print_global_state(qs);

    aqs::QCircuit fourier_gate = aqs::fourier_transform(qubits);
    qc << Gate(fourier_gate, 0);

    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "Output State (Fourier Basis):\n";
    aqs::print_global_state(qs);

    aqs::QCircuit inverse_fourier_gate = aqs::inverse_fourier_transform(qubits);
    qc.reset_circuit();
    qc << Gate(inverse_fourier_gate, 0);

    qc.generate_circuit();
    qs.simulate(qc);
    std::cout << "Output State after Inverse Fourier (Computational Basis):\n";
    print_global_state(qs);

    std::cout << "Output measurement after Fourier and Inverse Fourier: " << reverse_binary(qs.peek_measure_all(), qubits) << '\n';

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

void quantum_phase_estimation_example()
{
    using namespace aqs;
    std::cout << "---- Example: Phase Estimation Algorithm ----\n";

    int qcount = 7;
    aqs::QSimulator qs(qcount + 1);
    aqs::QCircuit qc(qcount + 1);

    qs.qubit(qcount) = aqs::QState::one();
    qs.generate_global_state();

    for (int i = 0; i < qcount; ++i)
        qc << H(i);

    float theta = 1.f/5.f;
    float angle = aqs::pi * 2.f * theta;
    
    std::cout << "\nEstimating a phase shift of " << theta << '\n';

    for (int i = 0; i < qcount; ++i)
        qc << CPhase(i, qcount, angle * static_cast<float>(1 << (qcount - 1 - i)));

    qc.generate_circuit();
    qs.simulate(qc);
    //std::cout << "Phase shifted state:\n";
    //aqs::print_global_state(qs);
    qs.generate_global_state();

    qc << Gate(aqs::inverse_fourier_transform(qcount), 0, "QFT†");

    qc.generate_circuit();
    qs.simulate(qc);
    //std::cout << "Phase estimation state output: ";
    //aqs::print_global_state(qs);

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
    std::cout << qc.representation() << "\n";

    std::cout << "-------------------------\n\n";
}

void quantum_shor_algorithm()
{
    using namespace aqs;
    std::cout << "---- Example: Shor's Algorithm ----\n";
    int output = 8;

    aqs::QCircuit h1(output / 2), h2(output - (output / 2));
    for (int i = 0; i < h1.qubit_count(); ++i)
        h1 << H(i);
    for (int i = 0; i < h2.qubit_count(); ++i)
        h2 << H(i);

    aqs::QSimulator qs(output + 4);
    aqs::QCircuit qc(output + 4);

    qs.qubit(qs.qubit_count() - 1) = aqs::QState::one();

    std::cout << "Adding Hadamard Gates...\n";
    qc << Gate(h1, 0);
    qc << Gate(h2, h1.qubit_count());


    // a^(2^j) mod 15 gate
    auto a2jmod15 = [](int a, int j) {
        aqs::QCircuit temp(4);
        if (a == 13 || a == 2)
        {
            temp << Swap(0, 1);
            temp << Swap(1, 2);
            temp << Swap(2, 3);
        }
        if (a == 7 || a == 8)
        {
            temp << Swap(2, 3);
            temp << Swap(1, 2);
            temp << Swap(0, 1);
        }
        if (a == 11 || a == 4)
        {
            temp << Swap(1, 3);
            temp << Swap(0, 2);
        }
        if (a == 7 || a == 11 || a == 13)
        {
            for (int j = 0; j < 4; ++j)
                temp << X(j);
        }

        for (int i = 0; i < j; ++i)
            temp << Gate(temp, 0);

        return temp;
    };

    int a = 13;
    int N = 15;
    for (int i = 0; i < output; ++i)
    {
        std::cout << "Adding Control_mod15 a^(2^" << i << ") Gate...\n";
        qc << ControlGate(a2jmod15(a, i), output - 1 - i, output, "a^(2^" + std::to_string(i) + ") mod 15");
    }

    std::cout << "Adding Fourier Transform Gate...\n";
    qc << Gate(aqs::inverse_fourier_transform(output), 0, "QFT†");

    std::cout << "Generating Global State...\n";
    qs.generate_global_state();

    std::cout << "Generating circuit...\n";
    qc.generate_circuit();

    std::cout << "Executing Simulation...\n";
    qs.simulate(qc);

    std::cout << "Profiling...\n";
    int reps = 1e4;
    auto profile = qs.profile_measure_all(reps);

    std::vector<float> qubit_probs(output);
    for (int i = 0; i < qubit_probs.size(); ++i)
        qubit_probs[i] = qs.qubit_probability_true(i);

    std::vector<std::pair<int, float>> state_probs(fast_pow2(output), {0 , 1.0f});
    for (int i = 0; i < state_probs.size(); ++i)
    {
        state_probs[i].first = i;
        for (int j = 0; j < output; ++j)
            state_probs[i].second *= (i & (1 << j)) ? qubit_probs[j] : (1.0f - qubit_probs[j]);
    }

    std::sort(state_probs.begin(), state_probs.end(), [](const std::pair<int, float>& a, const std::pair<int, float>& b){
        return b.second < a.second;
    });

    std::set<int> r_values;

    int top_count = 4;
    std::cout << "\nTop " << top_count << " Results:\n";
    for (int i = 0; i < top_count; ++i)
    {
        auto state = state_probs[i].first;
        auto prob = state_probs[i].second;
        auto base = fast_pow2(output);
        auto phase = static_cast<double>(state_probs[i].first) / static_cast<double>(base);
        auto frac = approximate_fraction(phase, 15);
        r_values.insert(static_cast<int>(frac.second));

        std::cout << binary_string(state, output) << " = " << std::setw(4) << state << " (Prob: " << std::setprecision(3) << prob
                  << std::setprecision(7) << ") - Phase: " << std::setw(4) << std::right << state << " / " << base << " = "
                  << std::setw(4) << phase << " = " << frac.first << " / " << frac.second << "\n";
    }


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

    //std::ofstream file;
    //file.open("out.txt");
    //file << aqs::gen_circuit_text_image(qc, qs);
    //file.close();

    std::cout << "-------------------------\n\n";
}

void quantum_constant_addition()
{
    using namespace aqs;
    std::cout << "---- Example: Constant Adder ----\n";
    int qubits = 5;
    int const_val = 8;
    int input_value = 1;

    aqs::QSimulator qs(qubits);
    aqs::QCircuit qc(qubits);
    
    qc << Gate(aqs::fourier_transform(qubits), 0, "QFT");

    aqs::QCircuit add_circuit(qubits);
    for (int i = 0; i < qubits; ++i)
        add_circuit << Phase(i, const_val * aqs::pi / static_cast<float>(1 << i));

    qc << Gate(add_circuit, 0);

    qc << Gate(aqs::inverse_fourier_transform(qubits), 0, "QFT†");
    for (int i = 0; i < qubits; ++i)
        qs.qubit(i) = (input_value & (1 << i)) ? aqs::QState::one() : aqs::QState::zero();

    qs.qubit(0) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::zero();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(4) = aqs::QState::zero();

    std::cout << "\nInput value: " << input_value << " (" << binary_string(input_value, qubits) << ") ; Constant Value: " << const_val << "\n";
    qs.generate_global_state();
    //std::cout << "Input state:\n";
    //aqs::print_global_state(qs);

    qc.generate_circuit();
    qs.simulate(qc);

    //std::cout << "Output state:\n";
    //aqs::print_global_state(qs);
    //aqs::print_profile(qs.profile_measure_all(1e4));

    std::cout << "Output value: " << binary_string(reverse_binary(qs.peek_measure_all(), qubits), qubits) << "\n";
    std::cout << "Result: " << input_value << " + " << const_val << " = " << reverse_binary(qs.peek_measure_all(), qubits) << "\n\n";

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "-------------------------\n\n";
}

void quantum_two_addition()
{
    using namespace aqs;
    std::cout << "---- Example: Two Number Adder ----\n";
    int input = 4;
    aqs::QSimulator qs(2 * input + 1);
    aqs::QCircuit qc(2 * input + 1);

    //Calculate |a> + |b>
    // Execute Fourier Transform on b input
    qc << Gate(aqs::fourier_transform(input + 1), input, "QFT");

    // Execute Fourier Transform on a input stored on b input
    for (int i = 0; i < input; ++i)
        qc << CPhase(i, input * 2, aqs::pi / (1 << (input - i)));

    for (int i = input; i >= 0; --i)
    {
        for (int j = 0; j < i; ++j)
            qc << CPhase(j, input + i - 1, aqs::pi / (1 << (i - j - 1)));
    }

    // Return back to computational basis
    qc << Gate(aqs::inverse_fourier_transform(input + 1), input, "QFT†");

    std::cout << "\nAddition of two numbers of " << input << " bits\n"
              << "First " << input << " qubits are number A and the next " << input << " qubits are number B\n"
              << "Algorithm computes C = A + B, and stores C in the last " << input + 1 << " qubits\n\n";

    for (int i = 0; i < fast_pow2(input * 2); ++i)
    {
        int valA = i % fast_pow2(input);
        int valB = i / fast_pow2(input);

        for (int i = 0; i < input; ++i)
            qs.qubit(i) = (valA & (1 << i)) ? aqs::QState::one() : aqs::QState::zero();
        for (int i = 0; i < input; ++i)
            qs.qubit(i + input) = (valB & (1 << i)) ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(2 * input) = aqs::QState::zero();

        qs.generate_global_state();

        qc.generate_circuit();

        qs.simulate(qc);

        int out_val = reverse_binary(qs.peek_measure_all() & ((1 << (input + 1)) - 1), input + 1);
        std::string out = binary_string(qs.peek_measure_all(), input + 1);
        std::reverse(out.begin(), out.end());
        
        std::cout << binary_string(valA, input) << " + " << binary_string(valB, input) << " = " << out
                  << " ( " << std::setw(2) << valA << " + " << std::setw(2) << valB << " = " << std::setw(2) << out_val << " )\n";
    }

    aqs::print_circuit_text_image(qc, qs);

    std::cout << "\n-------------------------\n\n";
}

void quantum_counting()
{
    std::cout << "---- Example: Quantum Grover Counting ----\n";

    using namespace aqs;
    int output_count = 6;
    int search_count = 2;

    aqs::QSimulator qs(output_count + search_count);
    aqs::QCircuit qc(output_count + search_count);

    for (int i = 0; i < output_count + search_count; ++i)
        qc << H(i);

    QCircuit oracle = grover_oracle(search_count, 0);
    for (int i = output_count - 1; i >= 0; --i)
        qc << ControlGate(aqs::grover_iteration(search_count, oracle, 1 << (output_count - 1 - i)), i, output_count);

    qc << Gate(aqs::inverse_fourier_transform(output_count), 0);
    qc.generate_circuit();
    qs.simulate(qc);

    std::vector<float> qubit_probs(output_count);
    for (int i = 0; i < qubit_probs.size(); ++i)
        qubit_probs[i] = qs.qubit_probability_true(i);

    std::vector<std::pair<uint32_t, float>> state_probs(fast_pow2(output_count), {0 , 1.0f});
    for (int i = 0; i < state_probs.size(); ++i)
    {
        state_probs[i].first = i;
        for (int j = 0; j < output_count; ++j)
            state_probs[i].second *= (i & (1 << j)) ? qubit_probs[j] : (1.0f - qubit_probs[j]);
    }

    std::sort(state_probs.begin(), state_probs.end(), [](const std::pair<uint32_t, float>& a, const std::pair<uint32_t, float>& b){
        return b.second < a.second;
    });

    auto likely_value = state_probs[0].first;
    auto phase = (float)likely_value * aqs::pi * 2.f / (float)(1 << output_count);

    auto search_space = 1 << search_count;
    auto solution_count = (1 << search_count) * std::cos(phase / 2.f) * std::cos(phase / 2.f);
    std::cout << "Likely output from phase estimation: " << likely_value << '\n';
    std::cout << "Estimated phase: " << phase << '\n';
    std::cout << "Estimated number of solutions: " << solution_count << '\n';

    auto error = (std::sqrt(2.f * solution_count * (search_space - solution_count)) + (float)search_space / (1 << output_count)) * 2.f / (1 << output_count);
    std::cout << "Error: " << error << "\n";

    for (int i = 0; i < output_count; ++i)
        std::cout << "i = " << i << " : " << qs.qubit_probability_true(i) << std::endl;

    std::cout << "\n-------------------------\n\n";
}

void temp()
{
    using namespace aqs;
    QCircuit qc(4);

    qc << Z{qc.qubit_count() - 1};

    for (int i = 0; i < qc.qubit_count(); ++i)
        qc << X(i);

    qc << Z{qc.qubit_count() - 1};

    qc << Gate(NControl_Gate(qc.qubit_count(), 0, qc.qubit_count() - 1, qc.qubit_count() - 1, Z::gate()), 0);

    for (int i = 0; i < qc.qubit_count(); ++i)
        qc << X(i);

    qc.generate_circuit();
    aqs::print_circuit_matrix(qc);

    QCircuit temp(1);
    temp << Y{0} << X{0} << Y{0} << Z{0} << X{0};
    temp << X{0} << Z{0} << Y{0} << X{0} << Y{0};
    temp.generate_circuit();
    temp.reset_circuit();
    temp << RotZ(0, aqs::pi * 2.f);
    temp.generate_circuit();
    aqs::print_circuit_matrix(temp);
}

void deutsch_jozsa()
{
    aqs::QCircuit qc(5);
    aqs::QSimulator qs(qc.qubit_count());

    qc << aqs::X{qc.qubit_count()};
    for (int i = 0; i < qc.qubit_count() - 1; ++i)
        qc << aqs::H(i);
}

void quantum_teleportation()
{
    std::cout << "---- Example: Quantum Teleportation ----\n";

    aqs::QCircuit qc(3);
    aqs::QSimulator qs(qc.qubit_count());

    float angle = -atan(0.75f);

    // Initial sender qubit state
    qs.qubit(0) = aqs::QState{0.6f, {0.0f, 0.8f}};
    qs.generate_global_state();

    std::cout << "Sender qubit:\n";
    aqs::print_state(qs.qubit(0));
    std::cout << '\n';

    // Create entangled pair
    qc << aqs::H{1} << aqs::CX{1, 2};
    qc << aqs::Barrier{};
    
    // Sender applies some operations
    qc << aqs::CX{0, 1} << aqs::H{0};

    // Sender measures their qubits
    qc.generate_circuit();
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

    qc0.generate_circuit();
    qs0.simulate(qc0);

    std::cout << "Profiling the qubit received:\n";
    aqs::print_profile(qs.profile_measure(2, 1e4));


    //This is to test if the state is the same as the sender's state
    qc0 << aqs::RotX{2, acos(0.6f) * 2};
    qc0.generate_circuit();
    qs.simulate(qc0);

    std::cout << "\nTesting qubit with a gate:\n";
    aqs::print_profile(qs.profile_measure(2, 1e4));
    std::cout << "Expected qubit measurement with gate operation: |0>\n";

    std::cout << "\n-------------------------\n\n";
}

void quantum_simulation()
{
    auto qc = aqs::QCircuit(2);
    auto qs = aqs::QSimulator(2);

    float theta, phi, lamda;
    auto ucirc = aqs::QCircuit(1);
    ucirc << aqs::Phase{0, lamda} << aqs::RotY{0, theta} << aqs::Phase{0, phi};
}

void quantum_hamiltonian()
{
    auto qc = aqs::QCircuit{4};
    auto qs = aqs::QSimulator{4};

    float timestep = 1e-3;

    auto hamiltonian = aqs::QCircuit{4};

    for (uint32_t i = 0; i < qc.qubit_count() - 1; ++i)
       hamiltonian << aqs::CX{i, qc.qubit_count() - 1};
    
    hamiltonian << aqs::RotZ{qc.qubit_count() - 1, aqs::pi - timestep};

    for (uint32_t i = 0; i < qc.qubit_count() - 1; ++i)
       hamiltonian << aqs::CX{qc.qubit_count() - 2 - i, qc.qubit_count() - 1};

    int steps = 1024;
    for (int i = 0; i < fast_log2(steps); ++i)
        hamiltonian << aqs::Gate(hamiltonian, 0, "Hamiltonian");

    qc << aqs::Gate{hamiltonian, 0};

    qc.generate_circuit();
    //aqs::print_circuit_matrix(qc);
    qs.simulate(qc);

    aqs::print_global_state(qs);
    aqs::print_circuit_text_image(qc, qs);

    aqs::QCircuit q{2};
    q << aqs::H{0} << aqs::CX{0, 1} << aqs::H{0};
    q.generate_circuit();
    aqs::QSimulator s{2};
    s.simulate(q);
    aqs::print_global_state(s);
    aqs::print_circuit_matrix(q);
    aqs::print_circuit_text_image(q, s);
}

void quantum_change_of_basis()
{
    aqs::QCircuit qc(2);
    qc  << aqs::X{0} << aqs::H{0} << aqs::H{1};
    qc.generate_circuit();
    aqs::QSimulator qs(2);
    qs.simulate(qc);
    aqs::print_global_state(qs);
    qs.set_basis(aqs::QSimulator::Basis::X);
    aqs::print_global_state(qs);
    qs.set_basis(aqs::QSimulator::Basis::Z);
    aqs::print_global_state(qs);
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << "\n";

    
    quantum_classic_xor_example();
    quantum_classic_or_example();
    quantum_2bit_adder_example();
    quantum_entanglement_example();
    quantum_superposition_example();
    quantum_grover_example();
    quantum_fourier_transform_example();

    quantum_phase_estimation_example();
    quantum_shor_algorithm();
    quantum_constant_addition();
    quantum_two_addition();

    quantum_teleportation();
    quantum_counting();


    quantum_change_of_basis();
   //quantum_hamiltonian();
    //temp();
}