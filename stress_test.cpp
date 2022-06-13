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
#include <vector>

#include "stress_test.h"

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

void X_gate_test(int qubit_count, int test_count);
void Y_gate_test(int qubit_count, int test_count);
void Z_gate_test(int qubit_count, int test_count);
void Phase_gate_test(int qubit_count, int test_count);
void swap_gate_test(int  qubit_count, int test_count);
void circuit_gate_test(int qubit_count, int gate_qcount, int test_count);
void profile_measure_all_test(int qubit_count, int test_count);
void profile_measure_test(int qubit_count, int test_count);

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << "\n";

    int repcount = 1000;
    int qcount = 10;
    int gateqcount = 4;
    
    //X_gate_test(qcount, repcount);

    //Y_gate_test(qcount, repcount);

    //Z_gate_test(qcount, repcount);

    Phase_gate_test(qcount, repcount);

    //circuit_gate_test(qcount, gateqcount, repcount);

    //profile_measure_all_test(qcount, repcount);

    //profile_measure_test(qcount, repcount);

    //swap_gate_test(qcount, repcount);

    return 0;
}

void X_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::X(0);
    };

    auto func2 = [&]() {
        qc << X_old(0);
    };

    profile(func1, test_count, "-- Test: New X gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, test_count, "-- Test: Old X gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Y_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Y(0);
    };

    auto func2 = [&]() {
        qc << Y_old(0);
    };

    profile(func1, test_count, "-- Test: New Y gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, test_count, "-- Test: Old Y gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Z_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Z(0);
    };

    auto func2 = [&]() {
        qc << Z_old(0);
    };

    profile(func1, test_count, "-- Test: New Z gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, test_count, "-- Test: Old Z gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void swap_gate_test(int  qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Swap(0, qubit_count - 1);
    };

    auto func2 = [&]() {
        qc << Swap_old(0, qubit_count - 1);
    };

    profile(func1, test_count, "-- Test: New Swap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, test_count, "-- Test: Old Swap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Phase_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    float angle = aqs::pi / 4.f;

    auto func1 = [&](){
        qc << aqs::Phase(0, angle);
    };

    auto func2 = [&]() {
        qc << Phase_old(0, angle);
    };

    profile(func1, test_count, "-- Test: New Phase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, test_count, "-- Test: Old Phase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
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

void profile_measure_all_test(int qubit_count, int test_count)
{
    aqs::QSimulator qs(qubit_count);

    auto func1 = [&]() {
        qs.profile_measure_all(test_count);
    };

    auto func2 = [&]() {
        std::vector<uint32_t> vals(qs.state_count());
        for (int i = 0; i < test_count; ++i)
            vals[qs.peek_measure_all()]++;
    };

    profile(func1, 1, "-- Test: New profile_all func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, 1, "-- Test: Old profile_all func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void profile_measure_test(int qubit_count, int test_count)
{
    aqs::QSimulator qs(qubit_count);
    uint32_t qubit_measured = 0;

    auto func1 = [&]() {
        qs.profile_measure(qubit_measured, test_count);
    };

    auto func2 = [&]() {
        std::array<uint32_t, 2> vals;
        uint32_t t = 0;
        for (int i = 0; i < test_count; ++i)
            t += qs.measure(qubit_measured);
        vals[0] = test_count - t;
        vals[1] = t;
    };

    profile(func1, 1, "-- Test: New profile func for circuit of " + std::to_string(qubit_count) + " qubits  --");

    profile(func2, 1, "-- Test: Old profile func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}