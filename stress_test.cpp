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


template<typename F>
void profile(F func, int test_count, const std::string& heading)
{
    using sys_clock = std::chrono::high_resolution_clock;
    af::sync();
    auto begin = sys_clock::now();
    for (int i = 0; i < test_count; ++i)
    {
        func();
    }
    af::sync();
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
void Hadamard_gate_test(int qubit_count, int test_count);
void swap_gate_test(int  qubit_count, int test_count);
void CX_gate_test(int qubit_count, int test_count);
void CY_gate_test(int qcount, int repcount);
void CZ_gate_test(int qcount, int repcount);
void CPhase_gate_test(int qcount, int repcount);
void CCNot_gate_test(int qcount, int repcount);
void Or_gate_test(int qcount, int repcount);
void CSwap_gate_test(int qcount, int repcount);
void CHadamard_gate_test(int qubit_count, int test_count);
void CircuitGate_test(int qubit_count, int gate_qcount, int test_count);
void ControlCircuitGate_test(int qubit_count, int gate_qcount, int test_count);
void circuit_gate_test(int qubit_count, int gate_qcount, int test_count);
void profile_measure_all_test(int qubit_count, int test_count);
void profile_measure_test(int qubit_count, int test_count);
void tensor_product_test(int size, int test_count);

af::array old_tensor(const af::array&, const af::array&);

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << "\n";

    int repcount = 1000;
    int qcount = 10;
    int gateqcount = 8;

    int size = 50;

    tensor_product_test(size, repcount);

    ControlCircuitGate_test(qcount, gateqcount, repcount);

    X_gate_test(qcount, repcount);

    Y_gate_test(qcount, repcount);

    Z_gate_test(qcount, repcount);

    Phase_gate_test(qcount, repcount);

    Hadamard_gate_test(qcount, repcount);

    CX_gate_test(qcount, repcount);

    CY_gate_test(qcount, repcount);

    CZ_gate_test(qcount, repcount);

    CPhase_gate_test(qcount, repcount);

    CCNot_gate_test(qcount, repcount);

    Or_gate_test(qcount, repcount);

    CSwap_gate_test(qcount, repcount);

    CHadamard_gate_test(qcount, repcount);

    CircuitGate_test(qcount, gateqcount, repcount);

    circuit_gate_test(qcount, gateqcount, repcount);

    profile_measure_all_test(qcount, repcount);

    profile_measure_test(qcount, repcount);

    swap_gate_test(qcount, repcount);

    return 0;
}

void X_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::X(0);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << X_old(0);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New X gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old X gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Y_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Y(0);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Y_old(0);
    };

    profile(func1, test_count, "-- Test: New Y gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Y gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Z_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Z(0);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Z_old(0);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New Z gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Z gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void swap_gate_test(int  qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Swap(0, qubit_count - 1);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Swap_old(0, qubit_count - 1);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New Swap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Swap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Phase_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    float angle = aqs::pi / 4.f;

    auto func1 = [&](){
        qc << aqs::Phase(0, angle);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Phase_old(0, angle);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New Phase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Phase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Hadamard_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Hadamard(0);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Hadamard_old(0);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New Hadamard gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Hadamard gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CX_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    uint32_t target = 0;
    uint32_t control = qubit_count - 1;

    auto func1 = [&](){
        qc << aqs::Control_X(control, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_X_old(control, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New ControlX gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old ControlX gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CY_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    uint32_t target = 0;
    uint32_t control = qubit_count - 1;

    auto func1 = [&](){
        qc << aqs::Control_Y(control, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_Y_old(control, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New ControlY gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old ControlY gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CZ_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    uint32_t target = 0;
    uint32_t control = qubit_count - 1;

    auto func1 = [&](){
        qc << aqs::Control_Z(control, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_Z_old(control, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New ControlZ gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old ControlZ gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CPhase_gate_test(int qubit_count, int test_count)
{
    aqs::QCircuit qc(qubit_count);

    uint32_t target = 0;
    uint32_t control = qubit_count - 1;
    float angle = aqs::pi / 4;

    auto func1 = [&](){
        qc << aqs::Control_Phase(control, target, angle);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_Phase_old(control, target, angle);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New ControlPhase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old ControlPhase gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CCNot_gate_test(int qubit_count, int test_count)
{
    int controlA = 0;
    int controlB = 1;
    int target = qubit_count - 1;

    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::CControl_Not(controlA, controlB, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << CControl_Not_old(controlA, controlB, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New CControl Not gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old CControl Not gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void Or_gate_test(int qubit_count, int test_count)
{
    int controlA = 0;
    int controlB = 1;
    int target = qubit_count - 1;

    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Or(controlA, controlB, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Or_old(controlA, controlB, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New Or gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old Or gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CSwap_gate_test(int qubit_count, int test_count)
{
    int controlA = 0;
    int controlB = 1;
    int target = qubit_count - 1;

    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Control_Swap(controlA, controlB, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_Swap_old(controlA, controlB, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New CSwap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old CSwap gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CHadamard_gate_test(int qubit_count, int test_count)
{
    int control = 0;
    int target = qubit_count - 1;

    aqs::QCircuit qc(qubit_count);

    auto func1 = [&](){
        qc << aqs::Control_Hadamard(control, target);
        qc.circuit().eval();
    };

    auto func2 = [&]() {
        qc << Control_Hadamard_old(control, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New CHadamard gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old CHadamard gate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void CircuitGate_test(int qubit_count, int gate_qcount, int test_count)
{
    int target = qubit_count - 1;

    aqs::QCircuit qc(qubit_count);
    aqs::QCircuit gate(gate_qcount);

    for (int i = 0; i < gate_qcount; ++i)
        gate << aqs::Hadamard(i);

    auto func1 = [&](){
        qc << aqs::CircuitGate(gate, 0);
        qc.circuit().eval();
    };

    auto func2 = [&](){
        qc << CircuitGate_old(gate, 0);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New CircuitGate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old CircuitGate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void ControlCircuitGate_test(int qubit_count, int gate_qcount, int test_count)
{
    int control = 0;
    int target = qubit_count - gate_qcount;

    aqs::QCircuit qc(qubit_count);
    aqs::QCircuit gate(gate_qcount);

    for (int i = 0; i < gate_qcount; ++i)
        gate << aqs::Hadamard(i);

    auto func1 = [&](){
        qc << aqs::ControlCircuitGate(gate, control, target);
        qc.circuit().eval();
    };

    auto func2 = [&](){
        qc << ControlCircuitGate_old(gate, control, target);
        qc.circuit().eval();
    };

    profile(func1, test_count, "-- Test: New ControlCircuitGate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
    qc.reset_circuit();
    profile(func2, test_count, "-- Test: Old ControlCircuitGate func for circuit of " + std::to_string(qubit_count) + " qubits  --");
}

void circuit_gate_test(int qubit_count, int gate_qcount, int test_count)
{
    aqs::QCircuit qc(qubit_count);
    aqs::QCircuit gate(gate_qcount);

    for (int i = 0; i < gate_qcount; ++i)
        gate << aqs::Hadamard(i);

    auto func = [&](){
        qc << aqs::CircuitGate(gate, 0);
        qc.circuit().eval();
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

void tensor_product_test(int size, int test_count)
{
    af::array a = af::randu(size, size, f32);
    af::array b = af::randu(size, size, f32);

    auto func1 = [&](){
        auto res = tensor_product(a, b);
        res.eval();
    };

    auto func2 = [&](){
        auto res = old_tensor(a, b);
        res.eval();
    };

    profile(func1, test_count, "-- Test: New tensor product for matrices of size " + std::to_string(size) + " --");

    profile(func2, test_count, "-- Test: Old tensor product for matrices of size " + std::to_string(size) + " --");
}

af::array old_tensor(const af::array& mat1, const af::array& mat2)
{
    af::dim4 dims1 = mat1.dims();
    af::dim4 dims2 = mat2.dims();

    if (mat1.type() != mat2.type() || dims1[2] != 1 || dims1[3] != 1 || dims2[2] != 1 || dims2[3] != 1)
        throw af::exception();

    af::array out = af::tile(mat2, dims1[0], dims1[1]);
    af::array resized_mat1 = af::resize(dims2[0], dims2[1], mat1, AF_INTERP_LOWER);

    return out *= resized_mat1;
}