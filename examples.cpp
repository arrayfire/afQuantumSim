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
    qc << CControl_Not(0, 2, 5);
    qc << CControl_Not(1, 3, 6);

    // Determine carry bits
    qc << Control_Not(0, 2);
    qc << Control_Not(1, 3);

    // Add carry bits from input
    qc << CControl_Not(3, 5, 6);

    // Add carry bits to output
    qc << Control_Not(2, 4);
    qc << Control_Not(3, 5);

    // Restore original input bits
    qc << Control_Not(0, 2);
    qc << Control_Not(1, 3);

    auto print_result = [](int val1, int val2, aqs::QSimulator& qs, aqs::QCircuit& qc) {
        qs.qubit(0) = val1 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(1) = val1 & 2 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(2) = val2 & 1 ? aqs::QState::one() : aqs::QState::zero();
        qs.qubit(3) = val2 & 2 ? aqs::QState::one() : aqs::QState::zero();

        qs.generate_global_state();
        std::cout << "Input: " << binary_string(qs.peek_measure_all(), 7) << '\n';

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
    qc << Hadamard(0);
    qc << Control_Not(0, 1);
    qs.simulate(qc);
    std::cout << "State of entangled q[0] and q[1]:\n";
    aqs::print_global_state(qs);
    //aqs::print_profile(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 3 qubits ***\n\n";
    //Entangle 3 qubits
    qs.generate_global_state();
    qc << Control_Not(0, 2);
    qs.simulate(qc);
    std::cout << "State of entangled q[0], q[1], and q[2]:\n";
    aqs::print_global_state(qs);
    //aqs::print_global_state(qs.profile_measure_all(reps));

    std::cout << "\n*** Entanglement of 4 qubits ***\n\n";
    //Entangle 4 qubits
    qs.generate_global_state();
    qc << Control_Not(0, 3);
    qs.simulate(qc);
    std::cout << "State of entangled q[0], q[1], q[2], and q[3]:\n";
    aqs::print_global_state(qs);
    //qs.profile_measure_all(reps);

    std::cout << "\n*** Entanglement of 2 pair of qubits ***\n\n";
    //Entangle 2 pair of qubits
    qs.generate_global_state();
    qc.reset_circuit();
    qc << Hadamard(0);
    qc << Hadamard(2);
    qc << Control_Not(0, 1);
    qc << Control_Not(2, 3);
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
    qc << Hadamard(0);
    qc << Hadamard(2);
    qs.simulate(qc);

    std::cout << "Global State with q[0] and q[2] in superposition using Hadamard:\n";
    print_global_state(qs);
    //print_profile_measure(qs.profile_measure_all(reps));

    // Put qubit 1 into a superposed state and qubit 3 into definite |1> state
    qs.qubit(1) = aqs::QState{{0.0f, 0.6f} , { 0.8f, 0.0f }};
    qs.qubit(3) = aqs::QState::one();
    qs.generate_global_state();
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
    qc << CircuitGate(grover_search_gate, 0);

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
    qc << CircuitGate(fourier_gate, 0);

    qs.simulate(qc);
    std::cout << "Output State (Fourier Basis):\n";
    aqs::print_global_state(qs);

    aqs::QCircuit inverse_fourier_gate = aqs::inverse_fourier_transform(qubits);
    qc.reset_circuit();
    qc << CircuitGate(inverse_fourier_gate, 0);

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
        qc << Hadamard(i);

    float theta = 1.f/5.f;
    float angle = aqs::pi * 2.f * theta;
    
    std::cout << "\nEstimating a phase shift of " << theta << '\n';

    for (int i = 0; i < qcount; ++i)
        qc << Control_Phase(i, qcount, angle * static_cast<float>(1 << (qcount - 1 - i)));

    qs.simulate(qc);
    //std::cout << "Phase shifted state:\n";
    //aqs::print_global_state(qs);
    qs.generate_global_state();

    qc << CircuitGate(aqs::inverse_fourier_transform(qcount), 0, "QFT†");
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
        h1 << Hadamard(i);
    for (int i = 0; i < h2.qubit_count(); ++i)
        h2 << Hadamard(i);

    aqs::QSimulator qs(output + 4);
    aqs::QCircuit qc(output + 4);

    qs.qubit(qs.qubit_count() - 1) = aqs::QState::one();

    std::cout << "Adding Hadamard Gates...\n";
    qc << CircuitGate(h1, 0);
    qc << CircuitGate(h2, h1.qubit_count());


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
            temp << CircuitGate(temp, 0);

        return temp;
    };

    int a = 13;
    int N = 15;
    for (int i = 0; i < output; ++i)
    {
        std::cout << "Adding Control_mod15 a^(2^" << i << ") Gate...\n";
        qc << ControlCircuitGate(a2jmod15(a, i), output - 1 - i, output, "a^(2^" + std::to_string(i) + ") mod 15");
    }

    std::cout << "Adding Fourier Transform Gate...\n";
    qc << CircuitGate(aqs::inverse_fourier_transform(output), 0, "QFT†");

    std::cout << "Generating Global State...\n";
    qs.generate_global_state();

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
    
    qc << CircuitGate(aqs::fourier_transform(qubits), 0, "QFT");

    aqs::QCircuit add_circuit(qubits);
    for (int i = 0; i < qubits; ++i)
        add_circuit << Phase(i, const_val * aqs::pi / static_cast<float>(1 << i));

    qc << CircuitGate(add_circuit, 0);

    qc << CircuitGate(aqs::inverse_fourier_transform(qubits), 0, "QFT†");
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
    qc << CircuitGate(aqs::fourier_transform(input + 1), input, "QFT");

    // Execute Fourier Transform on a input stored on b input
    for (int i = 0; i < input; ++i)
        qc << Control_Phase(i, input * 2, aqs::pi / (1 << (input - i)));

    for (int i = input; i >= 0; --i)
    {
        for (int j = 0; j < i; ++j)
            qc << Control_Phase(j, input + i - 1, aqs::pi / (1 << (i - j - 1)));
    }

    // Return back to computational basis
    qc << CircuitGate(aqs::inverse_fourier_transform(input + 1), input, "QFT†");

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
    int output_count = 3;
    int search_count = 2;

    aqs::QSimulator qs(output_count + search_count);
    aqs::QCircuit qc(output_count + search_count);

    for (int i = 0; i < output_count + search_count; ++i)
        qc << Hadamard(i);
    for (int i = 0; i < search_count; ++i)
        qc << ControlCircuitGate(aqs::grover_iteration(search_count, aqs::grover_oracle(search_count, 0), 1 << i), i, output_count);
    
    qc << CircuitGate(aqs::inverse_fourier_transform(output_count), 0);

    qs.simulate(qc);
    aqs::print_global_state(qs);

    //std::cout << qc.representation();
    //aqs::print_circuit_text_image(qc, qs);

    std::cout << "\n-------------------------\n\n";
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

    //quantum_counting();
}