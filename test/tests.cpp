/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#undef NDEBUG

#include "quantum.h"
#include "quantum_gates.h"
#include "quantum_visuals.h"

#include <iostream>

using namespace aqs;

bool fequal(float a, float b, float diff)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * diff);
};

bool cfequal(const af::cfloat& a, const af::cfloat& b, float diff)
{
    return fequal(a.real, b.real, diff) && fequal(a.imag, b.imag, diff);
}

bool qsequal(const aqs::QState& a, const aqs::QState& b, float diff)
{
    return cfequal(a[0], b[0], diff) && cfequal(a[1], b[1], diff);
}

bool arrayEquals(const af::array& lhs, const af::array& rhs, float abs_diff)
{
    return af::allTrue<bool>(af::abs(lhs - rhs) < abs_diff);
}
/**
 * @brief Tests all constructors of the QState class
 * 
 */
void test_qubit_constructor()
{
    std::cout << "Starting test_qubit_constructor..." << std::endl;

    aqs::QState q1;
    print_state(q1);
    assert((q1[0] == af::cfloat{ 1.f , 0.f } && q1[1] == af::cfloat{ 0.f , 0.f }));

    aqs::QState q2({0.6f, 0.f} , {0.f , 0.8f});
    print_state(q2);
    assert((q2[0] == af::cfloat{ 0.6f , 0.0f} && q2[1] == af::cfloat{ 0.f , 0.8f }));

    aqs::QState q3(q2);
    print_state(q3);
    assert((q3[0] == af::cfloat{ 0.6f , 0.f} && q3[1] == af::cfloat{0.f , 0.8f}));

    std::cout << "Finished test_qubit_constructor...\n" << std::endl;
}

/**
 * @brief Tests the assigment operator of the QState class
 * 
 */
void test_qubit_assigment()
{
    std::cout << "Starting test_qubit_assignment..." << std::endl;
    
    aqs::QState q1({0.6f, 0.f} , {0.f , 0.8f});
    aqs::QState q2;
    q2 = q1;

    print_state(q1);
    print_state(q2);

    assert((q2[0] == q1[0] && q2[1] == q1[1]));

    std::cout << "Finished test_qubit_assignment...\n" << std::endl;
}

/**
 * @brief Tests the normalization condition of the state stored by the QState class
 * 
 */
void test_qubit_normalization()
{
    std::cout << "Starting test_qubit_normalization..." << std::endl;

    aqs::QState q1({3.f, 0.f} , {4.f , 0.f});
    
    print_state(q1);
    assert((q1[0] == af::cfloat{ 0.6f , 0.0f} && q1[1] == af::cfloat{ 0.8f , 0.f }));

    aqs::QState q2({1.f, 2.f} , {3.f , 4.f});

    print_state(q2);
    af::cfloat norm_state2[2] = {
        {1.f , 2.f} ,
        {3.f , 4.f}
    };

    float mag = sqrt(norm_state2[0].real * norm_state2[0].real + norm_state2[0].imag * norm_state2[0].imag +
                     norm_state2[1].real * norm_state2[1].real + norm_state2[1].imag * norm_state2[1].imag);

    norm_state2[0] = norm_state2[0] / mag;
    norm_state2[1] = norm_state2[1] / mag;

    assert(q2[0] == norm_state2[0] && q2[1] == norm_state2[1]);
    
    std::cout << "Finished test_qubit_normalization...\n" << std::endl;
}

/**
 * @brief Tests the calculation of probability of each state by the QState class
 * 
 */
void test_qubit_probability()
{
    std::cout << "Starting test_qubit_probability..." << std::endl;

    float diff = 1e-4;

    aqs::QState q1({1.f , 0.f} , {0.f , 0.f});
    print_state(q1);
    assert(q1.probability_false() == 1.f);
    assert(q1.probability_true() == 0.f);

    aqs::QState q2({0.f , 1.f} , {0.f , 0.f});
    print_state(q2);
    assert(q2.probability_false() == 1.f);
    assert(q2.probability_true() == 0.f);

    aqs::QState q3({0.f , 0.f} , {1.f , 0.f});
    print_state(q3);
    assert(q3.probability_false() == 0.f);
    assert(q3.probability_true() == 1.f);

    aqs::QState q4({0.f , 0.f} , {0.f , 1.f});
    print_state(q4);
    assert(q4.probability_false() == 0.f);
    assert(q4.probability_true() == 1.f);

    aqs::QState q5({0.f , 1.f} , {0.f , 1.f});
    print_state(q5);
    assert(fequal(q5.probability_false(), 0.5f, diff));
    assert(fequal(q5.probability_true(), 0.5f, diff));

    aqs::QState q6({0.6f , 0.f} , {0.f , 0.8f});
    print_state(q6);
    assert(fequal(q6.probability_false(), 0.36f, diff));
    assert(fequal(q6.probability_true(), 0.64f, diff));

    std::cout << "Finished test_qubit_probability...\n" << std::endl;
}

/**
 * @brief Tests the randomization of the measurement (collapse) of states by the QState class
 * 
 */
void test_qubit_profile_measure()
{
    std::cout << "Starting test_qubit_profile_measure..." << std::endl;

    int reps = 1e4;

    //Produce 99.9% confidence interval
    auto range_true = [](const aqs::QState& q, int rep_count) {
        float std_dev = std::sqrt(q.probability_true() * q.probability_false());

        int expected = rep_count * q.probability_true();
        int diff = std::sqrt(rep_count) * 3.291f * std_dev;

        return std::pair<int, int>(expected - diff , expected + diff);
    };

    aqs::QState q1({1.f , 0.f} , {0.f , 0.f});
    auto q1prof = q1.profile_measure(reps);
    auto q1rng = range_true(q1, reps);
    assert((q1rng.first <= q1prof[1] && q1prof[1] <= q1rng.second));

    aqs::QState q2({0.f , 1.f} , {0.f , 0.f});
    auto q2prof = q2.profile_measure(reps);
    auto q2rng = range_true(q2, reps);
    assert((q2rng.first <= q2prof[1] && q2prof[1] <= q2rng.second));

    aqs::QState q3({0.f , 0.f} , {1.f , 0.f});
    auto q3prof = q3.profile_measure(reps);
    auto q3rng = range_true(q3, reps);
    assert((q3rng.first <= q3prof[1] && q3prof[1] <= q3rng.second));

    aqs::QState q4({0.f , 0.f} , {0.f , 1.f});
    auto q4prof = q4.profile_measure(reps);
    auto q4rng = range_true(q4, reps);
    assert((q4rng.first <= q4prof[1] && q4prof[1] <= q4rng.second));

    aqs::QState q5({0.f , 1.f} , {0.f , 1.f});
    auto q5prof = q5.profile_measure(reps);
    auto q5rng = range_true(q5, reps);
    assert((q5rng.first <= q5prof[1] && q5prof[1] <= q5rng.second));

    aqs::QState q6({0.6f , 0.f} , {0.f , 0.8f});
    auto q6prof = q6.profile_measure(reps);
    auto q6rng = range_true(q6, reps);
    assert((q6rng.first <= q6prof[1] && q6prof[1] <= q6rng.second));

    std::cout << "Finished test_qubit_profile_measure...\n" << std::endl;
}

/**
 * @brief Test qubits operation functions
 * 
 */
void test_qubit_operations()
{
    std::cout << "Starting test_qubit_operations..." << std::endl;

    float diff = 1e-4;

    const aqs::QState qref({0.5f, -0.5f}, {0.7f, 0.1f});
    print_state(qref);

    aqs::QState q1 = X_op(qref);
    print_state(q1);
    assert((q1[0] == af::cfloat{0.7f, 0.1f} && q1[1] == af::cfloat{0.5f, -0.5f}));

    aqs::QState q2 = Y_op(qref);
    print_state(q2);
    
    assert((q2[0] == af::cfloat{0.1f, -0.7f} && q2[1] == af::cfloat{0.5, 0.5f}));

    aqs::QState q3 = Z_op(qref);
    print_state(q3);
    
    assert((q3[0] == af::cfloat{0.5f, -0.5f} && q3[1] == af::cfloat{-0.7f, -0.1f}));

    aqs::QState q4 = Hadamard_op(qref);
    print_state(q4);

    aqs::QState q4h({1.2f,-0.4f},{-0.2f,-0.6f});

    assert((fequal(q4[0].real, q4h[0].real, diff) && fequal(q4[0].imag, q4h[0].imag, diff) &&
            fequal(q4[1].real, q4h[1].real, diff) && fequal(q4[1].imag, q4h[1].imag, diff)));

    std::cout << "Finished test_qubit_operations...\n" << std::endl;
}

/**
 * @brief Tests the different contructors of QSimulator
 * 
 */
void test_qsim_constructor()
{
    std::cout << "Starting test_qsim_constructor..." << std::endl;

    int count = 3;
    aqs::QState state{ {0.6f , 0.0f} , {0.0f , 0.8f} };
    aqs::QSimulator qs2(count, state);

    print_statevector(qs2);

    assert(qs2.qubit_count() == count);

    for (int i = 0; i < count; ++i)
        assert(qs2.qubit(i) == state);
    
    std::cout << "Finished test_qsim_constructor...\n" << std::endl;
}


/**
 * @brief Tests the collapsing of the statevector after measurement
 * 
 */
void test_qsim_measurement()
{
    std::cout << "Starting test_qsim_measurement..." << std::endl;
    int count = 3;

    aqs::QCircuit qc(count);
    aqs::QSimulator qs(count, QState::zero());

    qs.generate_statevector();
    print_statevector(qs);

    assert(qs.peek_measure_all() == 0);

    for (int i = 0; i < count; ++i)
        assert(qs.peek_measure(i) == false);

    qs.qubit(0) = QState::one();
    qs.qubit(count - 1) = QState::one();

    qs.generate_statevector();
    print_statevector(qs);

    assert(qs.peek_measure_all() == ((1 << (count - 1)) | 1));

    assert(qs.peek_measure(0) == true);
    assert(qs.peek_measure(count - 1) == true);
    for (int i = 1; i < count - 1; ++i)
        assert(qs.peek_measure(i) == false);

    qs.qubit(0) = aqs::QState({0.0f, 1.0f} , {-1.0f, 0.0f});
    qs.qubit(count - 1) = aqs::QState({0.6f, 0.0f}, {0.0f, 0.8f});
    qs.generate_statevector();

    print_statevector(qs);
    auto statevector = qs.measure_all();
    print_statevector(qs);

    assert((qs.state(statevector) == af::cfloat{1.0f, 0.0}));

    qs.qubit(0) = aqs::QState({0.0f, 1.0f} , {-1.0f, 0.0f});
    qs.qubit(count - 1) = aqs::QState({0.6f, 0.0f}, {0.0f, 0.8f});
    qs.generate_statevector();
    
    float diff = 1e-4;
    bool q0state = qs.measure(0);
    statevector = qs.peek_measure_all();
    assert((q0state == static_cast<bool>(statevector & (1 << (qs.qubit_count() - 1)))));

    print_statevector(qs);
    if (q0state)
    {
        assert(cfequal(qs.state(0b100), af::cfloat{-0.6f, 0.0f}, diff));
        assert(cfequal(qs.state(0b101), af::cfloat{0.0f, -0.8f}, diff));
    }
    else
    {
        assert(cfequal(qs.state(0b000), af::cfloat{0.0f, 0.6f}, diff));
        assert(cfequal(qs.state(0b001), af::cfloat{-0.8f, 0.0f}, diff));
    }

    std::cout << "Finished test_qsim_measurement...\n" << std::endl;

}


/**
 * @brief Tests the random generation of the measurements according to each state probability
 * 
 */
void test_qsim_profile_measure()
{
    std::cout << "Starting test_qsim_profile_measure..." << std::endl;

    int reps = 1e4;
    int qcount = 5;

    aqs::QCircuit qc(qcount);
    aqs::QSimulator qs(qcount);

    qs.qubit(3) = QState::one();
    qs.generate_statevector();

    // 99.9% confidence interval for statevector measurements
    auto state_range_true = [](const aqs::QSimulator& qs, int state, int rep_count) {
        float prob_true = qs.state_probability(state);
        float std_dev = std::sqrt(prob_true * (1.f - prob_true));

        int expected = rep_count * prob_true;
        int diff = std::sqrt(rep_count) * 3.291f * std_dev;

        return std::pair<int, int>(expected - diff , expected + diff);
    };

    // 99.9% confidence interval for qubit measurements
    auto qubit_range_true = [](const aqs::QState& q, int rep_count) {
        float std_dev = std::sqrt(q.probability_true() * q.probability_false());

        int expected = rep_count * q.probability_true();
        int diff = std::sqrt(rep_count) * 3.291f * std_dev;

        return std::pair<int, int>(expected - diff , expected + diff);
    };

    // Testing global measurements
    auto profile = qs.profile_measure_all(reps);
    for (int i = 0; i < profile.size(); ++i)
    {
        auto range = state_range_true(qs, i, reps);
        assert((range.first <= profile[i] && profile[i] <= range.second));
    }

    qs.qubit(0) = aqs::QState(1.0f, 1.0f);
    qs.qubit(1) = aqs::QState(-2.0f, 3.0f);
    qs.qubit(4) = aqs::QState({0.0f, 0.6f} , { 0.8f, 0.0f});
    qs.generate_statevector();
    profile = qs.profile_measure_all(reps);

    for (int i = 0; i < profile.size(); ++i)
    {
        auto range = state_range_true(qs, i, reps);
        assert((range.first <= profile[i] && profile[i] <= range.second));
    }

    // Testing individual qubits profile measure
    for (int i = 0; i < qcount; ++i)
    {
        auto prof = qs.profile_measure(i, reps);
        auto range = qubit_range_true(qs.qubit(i), reps);
        assert((range.first <= prof[1] && prof[1] <= range.second));
    }

    std::cout << "Finished test_qcircuit_profile_measure...\n" << std::endl;
}

/**
 * @brief Tests the functionality of the basic gates provided:
 * 
 * X, Y, Z, Hadamard, Phase, CNot, CCNot, CircuitGate, ControlCircuitGate
 * 
 */
void test_qsim_gates()
{
    std::cout << "Starting test_qsim_gates..." << std::endl;

    float diff = 1e-4;
    int qcount = 4;
    aqs::QSimulator qs(qcount);
    aqs::QCircuit qc(qcount);
    qs.qubit(3) = QState::one();

    // X gate test
    std::cout << "X gate\n";
    qs.generate_statevector();
    qc << aqs::X(2);
    qc << aqs::X(3);
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0010) == af::cfloat{1.0f, 0.f}));

    std::vector<aqs::X> xgates = { aqs::X{2} , aqs::X{3} };
    qc.clear();
    qs.generate_statevector();
    qc << xgates;
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0010) == af::cfloat{1.0f, 0.f}));

    // Y gate test
    std::cout << "Y gate\n";
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::one();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Y(2);
    qc << aqs::Y(3);
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0010) == af::cfloat{1.0f, 0.f}));

    std::vector<aqs::Y> ygates = { aqs::Y{2} , aqs::Y{3} };
    qc.clear();
    qs.generate_statevector();
    qc << ygates;
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0010) == af::cfloat{1.0f, 0.f}));

    // Z gate test
    std::cout << "Z gate\n";
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::one();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Z(2);
    qc << aqs::Z(3);
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0001) == af::cfloat{-1.0f, 0.f}));

    std::vector<aqs::Z> zgates = { aqs::Z{2} , aqs::Z{3} };
    qc.clear();
    qs.generate_statevector();
    qc << zgates;
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0001) == af::cfloat{-1.0f, 0.f}));

    // NOT gate test 
    std::cout << "Not gate\n";
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::one();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Not(2);
    qc << aqs::Not(3);
    qc.compile();
    qs.simulate(qc);
    assert((qs.state(0b0010) == af::cfloat{1.0f, 0.f}));

    // Hadamard gate test
    std::cout << "Hadamard gate\n";
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::one();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::H(2);
    qc << aqs::H(3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{-0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{-0.5f, 0.f}, diff));

    std::vector<aqs::H> hgates = { aqs::H{2} , aqs::H{3} };
    qc.clear();
    qs.generate_statevector();
    qc << hgates;
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{-0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{0.5f, 0.f}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{-0.5f, 0.f}, diff));

    // Phase gate test
    std::cout << "Phase gate\n";
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = QState::zero();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Phase(2, std::atan(0.75f));
    qc << aqs::Phase(3, std::atan(1.f));
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.6f, 0.f}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{0.f, 0.f}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{-0.48f, 0.64f}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{0.f, 0.f}, diff));

    std::vector<aqs::Phase> phasegates = { aqs::Phase{2 , std::atan(.75f) } , aqs::Phase{3 , std::atan(1.f)} };
    qc.clear();
    qs.generate_statevector();
    qc << phasegates;
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.6f, 0.f}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{0.f, 0.f}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{-0.48f, 0.64f}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{0.f, 0.f}, diff));

    // RotX gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::minus();
    qs.qubit(3) = aqs::QState::zero();

    qc << aqs::RotX(2, aqs::pi / 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.5, 0.5}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{0.0, 0.0}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{-0.5, -0.5}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{0.0, 0.0}, diff));

    // RotY gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::zero();

    qc << aqs::RotY(2, std::acos(0.6f) * 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{-0.8, 0.0}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{0.0, 0.0}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{0.6, 0.0}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{0.0, 0.0}, diff));

    // RotZ gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::minus();
    qs.qubit(3) = aqs::QState::zero();

    qc << aqs::RotZ(2, -aqs::pi / 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0000), af::cfloat{0.5, 0.5}, diff));
    assert(cfequal(qs.state(0b0001), af::cfloat{0.0, 0.0}, diff));
    assert(cfequal(qs.state(0b0010), af::cfloat{-0.5, 0.5}, diff));
    assert(cfequal(qs.state(0b0011), af::cfloat{0.0, 0.0}, diff));

    // XOR gate test
    std::cout << "Xor gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Xor(0, 2);
    qc << aqs::Xor(1, 3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.64f, 0.f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.f, 0.48f}, diff));

    // SWAP gate test
    std::cout << "Swap gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.0f, 0.8f} , {0.6f, 0.0f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Swap(0, 2);
    qc << aqs::Swap(1, 3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0001), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b1001), af::cfloat{-0.64f, 0.0f}, diff));
    assert(cfequal(qs.state(0b1101), af::cfloat{0.0f, 0.48f}, diff));

    // CNOT gate test
    std::cout << "CNOT gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << CNot(0, 2);
    qc << CNot(1, 3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.64f, 0.f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.f, 0.48f}, diff));

    // Control-X gate test
    std::cout << "CX gate\n";
    qs.qubit(0) = aqs::QState({1.0f, 0.0f} , {0.0f, 0.0f});
    qs.qubit(1) = aqs::QState({0.0f, 0.0f} , {1.0f, 0.0f});
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CX(0, 2);
    qc << aqs::CX(1, 3);
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.64f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.0f, 0.48f}, diff));

    //Group Control-X gate insertion test
    std::cout << "CX gate through group insertion\n";
    std::vector<aqs::CX> cxgates = { aqs::CX{0, 2} , aqs::CX{1, 3} };
    qs.generate_statevector();
    qc.clear();
    qc << cxgates;
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.64f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.0f, 0.48f}, diff));

    // Control-Y gate test
    std::cout << "CY gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CY(0, 2);
    qc << aqs::CY(1, 3);
    qc.compile();
    qs.simulate(qc);
   
    assert(cfequal(qs.state(0b0100), af::cfloat{0.48f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.0f, 0.36f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{0.0f, 0.64f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{-0.48f, 0.0f}, diff));

    // Control-Z gate test
    std::cout << "CZ gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CZ(0, 2);
    qc << aqs::CZ(1, 3);
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0100), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.0f, -0.48f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.64f, 0.0f}, diff));

    // Control-Phase gate test
    std::cout << "CPhase gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CPhase(0, 2, std::atan(0.75f));
    qc << aqs::CPhase(1, 3, std::atan(0.75f));
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.36f, 0.f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{-0.288f, 0.384f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{-0.512f, -0.384f}, diff));

    // CRotX gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::one();
    qs.qubit(2) = aqs::QState::minus();
    qs.qubit(3) = aqs::QState::minus();

    const float invsqrt8 = 0.3535533906f;
    qc << aqs::CRotX(0, 2, aqs::pi / 2.f);
    qc << aqs::CRotX(1, 3, aqs::pi / 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    print_statevector(qs);

    assert(cfequal(qs.state(0b0100), af::cfloat{invsqrt8, invsqrt8}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{-invsqrt8, -invsqrt8}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-invsqrt8, -invsqrt8}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{invsqrt8, invsqrt8}, diff));

    // CRotY gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::one();
    qs.qubit(2) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::one();

    qc << aqs::CRotY(0, 2, std::acos(0.6f) * 2.f);
    qc << aqs::CRotY(1, 3, std::acos(0.6f) * 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.0, 0.0}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.0, 0.0}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.8, 0.0}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.6, 0.0}, diff));

    // CRotZ gate test
    qc.clear();
    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::one();
    qs.qubit(2) = aqs::QState::minus();
    qs.qubit(3) = aqs::QState::minus();

    qc << aqs::CRotZ(0, 2, -aqs::pi / 2.f);
    qc << aqs::CRotZ(1, 3, -aqs::pi / 2.f);

    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{invsqrt8, invsqrt8}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{-invsqrt8, invsqrt8}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-invsqrt8, -invsqrt8}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{invsqrt8, -invsqrt8}, diff));

    // OR gate test
    std::cout << "OR gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::zero();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::Or(0, 1, 2);
    qc << aqs::Or(1, 0, 3);
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0111), af::cfloat{1.0f, 0.0f}, diff));

    qs.qubit(0) = QState::one();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b1111), af::cfloat{1.0f, 0.0f}, diff));
    
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::zero();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0000), af::cfloat{1.0f, 0.0f}, diff));

    // AND gate test
    std::cout << "And gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::zero();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::And(0, 1, 2);
    qc << aqs::And(1, 0, 3);
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0100), af::cfloat{1.0f, 0.0f}, diff));

    qs.qubit(0) = QState::one();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b1111), af::cfloat{1.0f, 0.0f}, diff));
    
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::zero();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0000), af::cfloat{1.0f, 0.0f}, diff));

    // CSWAP gate test
    std::cout << "CSwap gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = aqs::QState({0.6f, 0.0f} , {0.0f, 0.8f});
    qs.qubit(3) = aqs::QState({0.0f, 0.8f} , {0.6f, 0.0f});
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CSwap(0, 2, 3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{-0.64f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.0f, 0.48f}, diff));

    qs.generate_statevector();
    qc.clear();
    qc << aqs::CSwap(1, 2, 3);
    qc.compile();
    qs.simulate(qc);

    assert(cfequal(qs.state(0b0100), af::cfloat{0.0f, 0.48f}, diff));
    assert(cfequal(qs.state(0b0101), af::cfloat{-0.64f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0110), af::cfloat{0.36f, 0.0f}, diff));
    assert(cfequal(qs.state(0b0111), af::cfloat{0.0f, 0.48f}, diff));

    // CCNOT gate test
    std::cout << "CCNot gate\n";
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::one();
    qs.qubit(2) = QState::zero();
    qs.qubit(3) = QState::zero();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CCNot(0, 1, 2);
    qc << aqs::CCNot(1, 0, 3);
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0100), af::cfloat{1.0f, 0.0f}, diff));

    qs.qubit(0) = QState::one();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b1111), af::cfloat{1.0f, 0.0f}, diff));
    
    qs.qubit(0) = QState::zero();
    qs.qubit(1) = QState::zero();
    qs.generate_statevector();
    qc.compile();
    qs.simulate(qc);
    assert(cfequal(qs.state(0b0000), af::cfloat{1.0f, 0.0f}, diff));

    // Custom Control Gate
    std::cout << "Control Custom gate\n";
    aqs::QCircuit temp(2);
    temp << aqs::H(0);
    temp << aqs::CNot(0, 1);
    temp << aqs::Z(0);
    temp << aqs::Swap(0, 1);

    qc.clear();
    qc << ControlGate(temp, 0, 2);

    aqs::QCircuit ref(4);
    ref << aqs::CH(0, 2);
    ref << aqs::CCNot(0, 2, 3);
    ref << aqs::CZ(0, 2);
    ref << aqs::CSwap(0, 2, 3);

    qc.compile();
    ref.compile();
    assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));

    qc.clear();
    qc << ControlGate(temp, 2, 0);
    
    ref.clear();
    ref << aqs::CH(2, 0);
    ref << aqs::CCNot(0, 2, 1);
    ref << aqs::CZ(2, 0);
    ref << aqs::CSwap(2, 0, 1);
    
    qc.compile();
    ref.compile();
    assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));

    // Custom Gate
    std::cout << "Custom gate\n";
    qc.clear();
    temp.clear();
    temp << aqs::H(0);
    temp << aqs::CPhase(0, 1, 1.0f);
    temp << aqs::X(1);
    temp << aqs::Swap(0, 1);
    qc << Gate(temp, 1);

    ref.clear();
    ref << aqs::H(1);
    ref << aqs::CPhase(1, 2, 1.0f);
    ref << aqs::X(2);
    ref << aqs::Swap(1, 2);

    qc.compile();
    ref.compile();
    assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));

    // CHadamard Gate
    std::cout << "CHadamard gate\n";
    aqs::QSimulator reference(4);
    reference.qubit(0) = aqs::QState::zero();
    reference.qubit(3) = aqs::QState::zero();
    reference.qubit(1) = aqs::QState::zero();
    reference.qubit(2) = aqs::QState::one();
    reference.generate_statevector();

    qs.qubit(0) = aqs::QState::zero();
    qs.qubit(3) = aqs::QState::zero();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_statevector();
    qc.clear();
    qc << aqs::CH(0, 1);
    qc << aqs::CH(3, 2);

    qc.compile();
    qs.simulate(qc);
    
    assert(af::allTrue<bool>(reference.statevector() == qs.statevector()));

    float rsqrt2 = std::sqrt(1.0f / 2.0f);
    reference.qubit(0) = aqs::QState::one();
    reference.qubit(3) = aqs::QState::one();
    reference.qubit(1) = aqs::QState{ rsqrt2 , rsqrt2 };
    reference.qubit(2) = aqs::QState{ rsqrt2 , -rsqrt2 };
    reference.generate_statevector();

    qs.qubit(0) = aqs::QState::one();
    qs.qubit(3) = aqs::QState::one();
    qs.qubit(1) = aqs::QState::zero();
    qs.qubit(2) = aqs::QState::one();
    qs.generate_statevector();

    qc.compile();
    qs.simulate(qc);
    
    assert(af::allTrue<bool>(reference.statevector() == qs.statevector()));

    std::cout << "Finished test_qsim_gates...\n" << std::endl;
}

/**
 * @brief Tests the special gates provided by quantum_gates.cpp
 * 
 */
void test_special_gates()
{
    std::cout << "Starting test_special_gates..\n";

    {
        QCircuit qc(6);
        QSimulator qs(6);

        uint32_t target = 3;

        qc << aqs::Gate(aqs::NControl_Gate(6, { 0, 2, 4, 5 }, target, X::gate()), 0);

        for (int i = 0; i < qs.qubit_count(); ++i)
            qs.qubit(i) = aqs::QState::one();

        qc.compile();
        qs.generate_statevector();
        qs.simulate(qc);
        assert(qs.peek_measure_all() == 0b111011);

        qs.qubit(target) = aqs::QState::zero();
        qs.generate_statevector();
        qs.simulate(qc);
        assert(qs.peek_measure_all() == 0b111111);

        qs.qubit(1) = aqs::QState::zero();
        qs.qubit(target) = aqs::QState::zero();
        qs.generate_statevector();
        qs.simulate(qc);
        assert(qs.peek_measure_all() == 0b101111);

        qs.qubit(1) = aqs::QState::zero();
        qs.qubit(target) = aqs::QState::one();
        qs.generate_statevector();
        qs.simulate(qc);
        assert(qs.peek_measure_all() == 0b101011);
    }

    {
        QCircuit qc(4);
        qc << Gate(NControl_Gate(4, {3}, 0, X::gate()), 0);

        QCircuit ref(4);
        ref << CX(3, 0);

        qc.compile();
        ref.compile();

        assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));

        qc.clear();
        ref.clear();

        qc << Gate(NControl_Gate(4, {0, 3}, 2, X::gate()), 0);
        ref << CCNot(0, 3, 2);

        qc.compile();
        ref.compile();

        assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));
    }

    {
        QCircuit qc(4);
        qc << Gate(Control_Group_Gate(4, 1, {0, 2, 3}, X::gate()), 0);
        
        QCircuit ref(4);
        ref << CX{1, 0} << CX{1, 2} << CX{1, 3};

        qc.compile();
        ref.compile();

        assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));
    }

    {
        QCircuit qc = Group_Gate(4, {0, 2, 3}, X::gate(), true);

        QCircuit ref{ 4 };
        ref << X{ 0 } << X{ 2 } << X{ 3 };
        ref.compile();

        assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));
    }

    {
        QCircuit circuit{ 4 };
        circuit << H{0} << CX{0 , 2} << H{1} << Z{3};
        circuit.compile();

        QCircuit qc = Rewire_Gate(4, { 2 , 1 , 3 , 0 }, circuit, true);

        QCircuit ref{ 4 };
        ref << H{2} << CX{ 2 , 3 } << H{1} << Z{0};
        ref.compile();

        assert(af::allTrue<bool>(qc.circuit() == ref.circuit()));
    }

    {
        QCircuit ref{ 4 };
        ref 
        << CRotY{ 3 , 0,  aqs::pi / 2.f } 
        << X{ 0 } 
        << CRotX{ 1 , 3, -aqs::pi / 4.f } 
        << CPhase{ 2 , 1,  aqs::pi / 3.f} 
        << Swap{ 1 , 2 };
        ref.compile();

        QCircuit circuit{ 4 };
        circuit 
        << Swap{ 1 , 2 } 
        << CPhase{ 2 , 1 , -aqs::pi / 3.f } 
        << CRotX{ 1 , 3 , aqs::pi / 4.f } 
        << X{ 0 } 
        << CRotY{ 3 , 0 , -aqs::pi / 2.f };
        circuit.compile();

        QCircuit qc = Adjoint_Gate(circuit);

        assert(arrayEquals(qc.circuit(), ref.circuit(), 1e-5));
    }

    std::cout << "Finished test_special_gates...\n" << std::endl;
}

/**
 * @brief Tests the calculation of probabilities by QSimulator
 * 
 */
void test_qsim_probability()
{
    std::cout << "Starting test_qsim_probability..." << std::endl;
    int qcount = 4;
    float diff = 1e-4;
    aqs::QSimulator qs(4);

    qs.generate_statevector();
    for (int i = 0; i < qcount; ++i)
    {
        assert(qs.qubit_probability_true(i) == 0.f);
        assert(qs.qubit_probability_false(i) == 1.f);
    }

    qs.qubit(qcount - 1) = aqs::QState{0.f , 1.f};
    qs.generate_statevector();

    for (int i = 0; i < qcount - 1; ++i)
    {
        assert(qs.qubit_probability_true(i) == 0.f);
        assert(qs.qubit_probability_false(i) == 1.f);
    }

    assert(qs.qubit_probability_true(qcount - 1) == 1.f);
    assert(qs.qubit_probability_false(qcount - 1) == 0.f);

    qs.qubit(qcount - 1) = aqs::QState{1.f , 1.f};
    qs.generate_statevector();

    for (int i = 0; i < qcount - 1; ++i)
    {
        assert(qs.qubit_probability_true(i) == 0.f);
        assert(qs.qubit_probability_false(i) == 1.f);
    }

    assert(fequal(qs.qubit_probability_true(qcount - 1), 0.5f, diff));
    assert(fequal(qs.qubit_probability_false(qcount - 1), 0.5f, diff));

    qs.qubit(0) = aqs::QState{{0.6f , 0.0f} , {0.f , 0.8f}};
    qs.generate_statevector();

    for (int i = 1; i < qcount - 1; ++i)
    {
        assert(qs.qubit_probability_true(i) == 0.f);
        assert(qs.qubit_probability_false(i) == 1.f);
    }

    assert(fequal(qs.qubit_probability_true(0), 0.64f, diff));
    assert(fequal(qs.qubit_probability_false(0), 0.36f, diff));
    assert(fequal(qs.qubit_probability_true(qcount - 1), 0.5f, diff));
    assert(fequal(qs.qubit_probability_false(qcount - 1), 0.5f, diff));

    std::cout << "Finished test_qsim_probability...\n" << std::endl;
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);

    test_qubit_constructor();
    test_qubit_assigment();
    test_qubit_normalization();
    test_qubit_probability();
    test_qubit_profile_measure();
    test_qubit_operations();

    test_qsim_constructor();
    test_qsim_measurement();
    test_qsim_profile_measure();
    test_qsim_gates();
    test_qsim_probability();
    test_special_gates();

    return 0;
}