/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include <iostream>
#include <chrono>
#include <string>
#include <vector>

#include "quantum.h"

using sys_clock = std::chrono::high_resolution_clock;

template<typename F>
void profile(F func, int test_count, const std::string& heading)
{
    auto begin = sys_clock::now();
    for (int i = 0; i < test_count; ++i)
        func();
    auto end = sys_clock::now();

    auto diff = end - begin;
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    auto avg_time = std::chrono::duration_cast<std::chrono::microseconds>(diff / test_count).count();
    std::cout << heading << "\ntest count: " << test_count <<  " ; total time: " << total_time <<  " ms ; avg time: " << avg_time << " us\n";
}

void one_qubit_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func = [&](){
        qc << aqs::Hadamard(0);
    };

    profile(func, test_count, "-- Test: " + std::to_string(qubit_count) + " qubits  --");
}

void circuit_gate_test(int qubit_count, int gate_qcount, int test_count)
{
    aqs::QCircuit qc(qubit_count);
    aqs::QCircuit gate(gate_qcount);

    for (int i = 0; i < gate_qcount; ++i)
        gate << aqs::Hadamard(i);

    auto func = [&](){
        qc << aqs::CircuitGate(gate, 0);
    };

    profile(func, test_count, "-- Test: Adding a Gate of " + std::to_string(gate_qcount) + " qubits to a circuit of " +
            std::to_string(qubit_count) + " qubits  --");
}

void profile_measure_test(int qubit_count, int test_count)
{
    aqs::QSimulator qs(qubit_count);

    auto func1 = [&]() {
        qs.profile_measure_all(test_count);
    };

    auto func2 = [&]() {
        std::vector<int> vals(qs.state_count());
        for (int i = 0; i < test_count; ++i)
            vals[qs.peek_measure_all()]++;
    };

    profile(func1, 1, "-- Test: New profile for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, 1, "-- Test: Old profile for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << "\n";

    int repcount = 10000;
    int qcount = 8;
    int gateqcount = 4;
    
    //one_qubit_gate_test(qcount, repcount);

    //circuit_gate_test(qcount, gateqcount, repcount);

    profile_measure_test(qcount, repcount);

    return 0;
}