/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#pragma once

#include <arrayfire.h>

#include <array>
#include <cassert>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "utils.h"

/**
 * @brief ArrayFire Quantum Simulator
 * 
 */
namespace aqs
{

static const float pi = 3.14159265358979323846f;
static constexpr uint32_t max_qubit_count = 30;

class QGate;
class QState;
class QCircuit;
class QSimulator;

/**
 * @brief Class managing the behaviour of one qubit
 * 
 */
class QState
{
public:
    QState() = default;
    QState(const QState&) = default;
    QState(QState&&) = default;

    ~QState() = default;

    QState& operator=(const QState&) = default;
    QState& operator=(QState&&) = default;

    /**
     * @brief Construct a new QState object given the coefficients of the state
     * 
     * @details It normalizes the state passed
     * 
     * @param zeroState coefficient of the |0> state
     * @param oneState coefficient of the |1> state
     */
    QState(const std::complex<float>& zeroState, const std::complex<float>& oneState);

    /**
     * @brief Construct a new QState object given the coefficients of the state
     * 
     * @details It normalizes the state passed
     * 
     * @param states array with the complex number state |0> in the [0] entry and |1> in the [1] entry
     */
    QState(const std::array<std::complex<float>, 2>& states);

    /**
     * @brief Create a qubit state object given its polar coordinates
     * 
     * @param polar_angle angle in the xy plane
     * @param azimuthal_angle angle with respect to the xy plane
     * @return std::array<std::complex<float>, 2> 
     */
    static std::array<std::complex<float>, 2> create_state(float polar_angle, float azimuthal_angle)
    {
        return { std::cos(polar_angle / 2.f) ,
                 std::complex<float>{ std::cos(azimuthal_angle),
                 std::sin(azimuthal_angle)} * std::sin(polar_angle / 2.f) };
    }
    
    /**
     * @brief Sets the state of the object to given the coefficients of the state
     * 
     * @details It normalizes the states passed
     * 
     * @param zero_state coefficient of the |0> state
     * @param one_state coefficient of the |1> state
     * @return QState& reference to current QState object
     */
    QState& set(const std::complex<float>& zero_state, const std::complex<float>& one_state);

    /**
     * @brief Executes a measurement without collapsing the state
     * 
     * @return true for |1> state measurement
     * @return false for |0> state measurement
     */
    bool peek_measure() const;

    /**
     * @brief Executes a measurement and collapses the state to it
     * 
     * @return true for |1> state measurement
     * @return false for |0> state measurement
     */
    bool measure();

    /**
     * @brief Executes the given number of measurements on the given state
     *        and returns the results of the number of measurements of each state
     * 
     * @param rep_count number of measurements to be done
     * @return std::array<int, 2> array containing the results
     */
    std::array<uint32_t, 2> profile_measure(uint32_t rep_count) const;

    /**
     * @brief Returns the probability of finding the state in the |1> state
     * 
     * @return float probability |1> measurement
     */
    float probability_true() const noexcept { return state_[1].real * state_[1].real + state_[1].imag * state_[1].imag; }

    /**
     * @brief Returns the probability of finding the state in the |0> state
     * 
     * @return float probability of |0> measurement
     */
    float probability_false() const noexcept { return 1.f - probability_true(); }

    /**
     * @brief Returns the complex number of the state selected
     * 
     * @param index 0 or 1 depending on the state to obtain
     * @return const af::cfloat& 
     */
    const af::cfloat& operator[](bool index) const noexcept { return state_[static_cast<int>(index)]; }

    bool operator!=(const aqs::QState& other) const noexcept { return state_[0] != other.state_[0] || state_[1] != other.state_[1]; }
    bool operator==(const aqs::QState& other) const noexcept { return !(*this != other); }

    /**
     * @brief Returns a pointer to the array of stored |0> and |1> states
     * 
     * @return af::cfloat* 
     */
    af::cfloat* data() noexcept { return state_; }

    /**
     * @brief Returns a const pointer to the array of stored |0> and |1> states
     * 
     * @return af::cfloat* 
     */
    const af::cfloat* data() const noexcept { return state_; }

    /**
     * @brief Returns an af::array of the stored |0> and |1> states
     * 
     * @return af::array with dims [2,1,1,1]
     */
    af::array to_array() const { return af::array(2, state_); }

    /**
     * @brief Returns a qubit in zero |0> state
     * 
     * @return const QState& 
     */
    static const QState& zero() { const static QState zero_state{ 1.0f, 0.0f }; return zero_state; };

    /**
     * @brief Returns a qubit in one |1> state
     * 
     * @return const QState& 
     */
    static const QState& one () { const static QState  one_state{ 0.0f, 1.0f }; return  one_state; };

    /**
     * @brief Returns a qubit in plus |+> state (1/sqrt(2) |0> + 1/sqrt(2) |1>)
     * 
     * @return const QState& 
     */
    static const QState& plus() { const static QState plus_state{ 0.70710678118f , 0.70710678118f }; return plus_state; }

    /**
     * @brief Returns a qubit in minus |-> state (1/sqrt(2) |0> - 1/sqrt(2) |1>)
     * 
     * @return const QState& 
     */
    static const QState& minus() { const static QState minus_state{ 0.70710678118f , -0.70710678118f }; return minus_state; }

private:
    af::cfloat state_[2] { { 1.0f , 0.0f } , { 0.0f , 0.0f } };

    /**
     * @brief Normalizes the internal states of the qubit
     * 
     */
    void force_normalize();
};

class QCircuit
{
public:
    QCircuit() = delete;
    QCircuit(const QCircuit&) = default;
    QCircuit(QCircuit&&) = default;

    ~QCircuit() = default;

    QCircuit& operator=(const QCircuit&) = default;
    QCircuit& operator=(QCircuit&&) = default;

    /**
     * @brief Construct a new QCircuit object given the number of qubits
     *        and initial state for all the qubits
     * 
     * @param qubit_count number of qubits for the simulation
     */
    QCircuit(uint32_t qubit_count);

    // Global Operations
    void Global_Measure();

    void Measure(uint32_t qubit);

    template<typename T>
    friend QCircuit& operator<<(QCircuit& qc, const T& gate);

    template<typename T>
    friend QCircuit& operator<<(QCircuit& qc, const std::vector<T>& gates);

    uint32_t qubit_count() const noexcept { return qubits_; }
    uint32_t state_count() const noexcept { return fast_pow2(qubits_); }

    af::array& circuit() noexcept { return circuit_; }
    const af::array& circuit() const noexcept { return circuit_; }

    auto& gate_list() noexcept { return gate_list_; }
    const auto& gate_list() const noexcept { return gate_list_; }

    friend class QSimulator;

    const std::string& representation() const noexcept { return representation_; }

    void generate_circuit();

    /**
     * @brief Cleans the circuit removing all gates
     * 
     */
    void reset_circuit();

private:
    std::vector<std::shared_ptr<QGate>> gate_list_;
    af::array circuit_;
    std::string representation_;
    uint32_t qubits_ = 0;
    std::size_t cached_index_ = 0;
};

class QNoise
{
};

class QSimulator
{
public:
    /**
     * @brief Construct a new QSimulator object
     * 
     * @param qubit_count number of qubits of the simulator
     * @param initial_state initial state for all the qubits
     * @param noise_generator the noise generator to be used
     */
    explicit QSimulator(uint32_t qubit_count,
                        const QState& initial_state   = aqs::QState::zero(),
                        const QNoise& noise_generator = QNoise{});

    /**
     * @brief Calculates the global state vector from all the individual states of the qubits
     * 
     */
    void generate_global_state();
    
    /**
     * @brief Updates the global state of the simulation using the computation from the circuit
     * 
     * @param circuit circuit to simulate
     * 
     */
    void simulate(const QCircuit& circuit);

    /**
     * @brief Measures the given qubit without collapsing the state
     * 
     * @param qubit index of the qubit
     * 
     * @return true for |1> state measurement
     * @return false for |0> state measurement
     */
    bool peek_measure(uint32_t qubit) const;

    /**
     * @brief Measures the given qubit and collapses its state (and updates the global states) to it
     * 
     * @details In the global state, all states that contain the given measurement are updated according to the previous probability
     *          using conditional probability, while the other states are set to 0
     * 
     * @param qubit index of the qubit
     * 
     * @return true for |1> state measurement
     * @return false for |0> state measurement
     */
    bool measure(uint32_t qubit);

    /**
     * @brief Returns a bit representation of the global state measured after measuring all bits.
     *        It collapses the global state to that state
     * 
     * @return uint32_t measurement
     */
    uint32_t measure_all();

    /**
     * @brief Returns a bit representation of the global state measured after measuring all bits.
     *        It does not collpase the global state
     * 
     * @return uint32_t measurement
     */
    uint32_t peek_measure_all() const;

    /**
     * @brief Profiles the measurement of the given state of a qubit from the stored global state using the given number of tests
     * 
     * @param qubit qubit to measure
     * @param rep_count number of measurements to be done for profiling
     * @return std::array<int, 2> [0] = # of measurements for state |0>; [1] = # of measurements for state |1>
     */
    std::array<uint32_t, 2> profile_measure(uint32_t qubit, uint32_t rep_count) const;

    /**
     * @brief Profiles the measurement of the given output state from the stored global state using the given number of tests
     * 
     * @param rep_count number of measurements to be done for profiling
     * @return std::vector<int> vector with the measurements [int(xxxx...)] = #number of measurements of state |xxxxx> where xxxxx is the binary output state
     */
    std::vector<uint32_t> profile_measure_all(uint32_t rep_count) const;

    QState& qubit(uint32_t index) noexcept { assert(index < qubit_count()); return states_[index]; };
    const QState& qubit(uint32_t index) const noexcept { assert(index < qubit_count()); return states_[index]; };

    /**
     * @brief Returns the complex number of the state in the global state vector
     * 
     * @param state int representation of the binary unique state
     * @return af::cfloat 
     */
    af::cfloat state(uint32_t state) const noexcept { assert(state < state_count()); return global_state_(state).scalar<af::cfloat>(); }

    /**
     * @brief Returns the probability to measure the given qubit in the |1> state
     * 
     * @param qubit index of the qubit
     * @return float probability of |1> state
     */
    float qubit_probability_true(uint32_t qubit) const;

    /**
     * @brief Returns the probability to measure the given qubit in the |0> state
     * 
     * @param qubit index of the qubit
     * @return float probability of |0> state
     */
    float qubit_probability_false(uint32_t qubit) const { return 1.f - qubit_probability_true(qubit); }

    /**
     * @brief Returns the probability to measure the given state from teh global state
     * 
     * @param state int representation of the binary unique state
     * @return float probability of |xxxxxx> where xxxxx is the state
     */
    float state_probability(uint32_t state) const;

    /**
     * @brief Returns the number of qubits of the simulator
     * 
     * @return int number of qubits
     */
    uint32_t qubit_count() const noexcept { return qubits_; }

    /**
     * @brief Returns the number of unique states that can be measured from the global state
     * 
     * @return int number of states
     */
    uint32_t state_count() const noexcept { return fast_pow2(qubits_);}

    /**
     * @brief Returns a const af::array& to the internal global state stored by the simulator
     * 
     * @return const af::array& 
     */
    const af::array& global_state() const noexcept { return global_state_; }

private:
    std::vector<QState> states_;
    af::array global_state_;
    QNoise noise_;
    uint32_t qubits_;
};

class QGate
{
public:
    virtual QCircuit& operator()(QCircuit& qc) const = 0;
    virtual std::string to_string() const = 0;
};

template<typename T>
QCircuit& operator<<(QCircuit& qc, const T& gate)
{
    static_assert(std::is_base_of<QGate, T>::value, "Gate must inherit from QGate class");
    qc.representation_.append(gate.to_string());
    //return gate(qc);
    qc.gate_list().push_back(std::make_shared<T>(gate));
    return qc;
}

template<typename T>
QCircuit& operator<<(QCircuit& qc, const std::vector<T>& gates)
{
    static_assert(std::is_base_of<QGate, T>::value, "Gate must inherit from QGate class");
    for (const auto& gate : gates)
    {
        qc.representation_.append(gate.to_string());
        //gate(qc);
        qc.gate_list().push_back(std::make_shared<T>(gate));
    }
    return qc;
}

class Barrier : public QGate
{
public:
    Barrier(bool visible_ = true) noexcept
        :visible{visible_} {}
    QCircuit& operator()(QCircuit& qc) const override { return qc; }
    std::string to_string() const override { return visible ? "B;" : "P;"; }
    bool visible = true;
};

class X : public QGate
{
public:
    X(uint32_t target_qubit_) : target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit;
};

using Not = X;

class Y : public QGate
{
public:
    Y(uint32_t target_qubit_) : target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit;
};

class Z : public QGate
{
public:
    Z(uint32_t target_qubit_) : target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit;
};

class Hadamard : public QGate
{
public:
    Hadamard(uint32_t target_qubit_) : target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit;
};

class Phase : public QGate
{
public:
    Phase(uint32_t target_qubit_, float angle_) : target_qubit{target_qubit_}, angle{angle_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit;
    float angle;
};

class Swap : public QGate
{
public:
    Swap(uint32_t target_qubit_A_, uint32_t target_qubit_B_) : target_qubit_A{target_qubit_A_},
                                                     target_qubit_B{target_qubit_B_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t target_qubit_A;
    uint32_t target_qubit_B;
};

class Control_X : public QGate
{
public:
    Control_X(uint32_t control_qubit_, uint32_t target_qubit_) : control_qubit{control_qubit_},
                                                       target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit;
};

class Control_Y : public QGate
{
public:
    Control_Y(uint32_t control_qubit_, uint32_t target_qubit_) : control_qubit{control_qubit_},
                                                       target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit;
};

class Control_Z : public QGate
{
public:
    Control_Z(uint32_t control_qubit_, uint32_t target_qubit_) : control_qubit{control_qubit_},
                                                       target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit;
};

class Control_Hadamard : public QGate
{
public:
    Control_Hadamard(uint32_t control_qubit_, uint32_t target_qubit_) : control_qubit{control_qubit_},
                                                              target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit;
};

using Control_Not = Control_X;
using Xor = Control_X;

class Control_Phase : public QGate
{
public:
    Control_Phase(uint32_t control_qubit_, uint32_t target_qubit_, float angle_) : control_qubit{control_qubit_},
                                                                         target_qubit{target_qubit_},
                                                                         angle{angle_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit;
    float angle;
};

class Control_Swap : public QGate
{
public:
    Control_Swap(uint32_t control_qubit_, uint32_t target_qubit_A_, uint32_t target_qubit_B_)
        : control_qubit{control_qubit_}, target_qubit_A{target_qubit_A_}, target_qubit_B{target_qubit_B_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit;
    uint32_t target_qubit_A;
    uint32_t target_qubit_B;
};

class CControl_Not : public QGate
{
public:
    CControl_Not(uint32_t control_qubit_A_, uint32_t control_qubit_B_, uint32_t target_qubit_)
        : control_qubit_A{control_qubit_A_}, control_qubit_B{control_qubit_B_}, target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit_A;
    uint32_t control_qubit_B;
    uint32_t target_qubit;
};

using And = CControl_Not;

class Or : public QGate
{
public:
    Or(uint32_t control_qubit_A_, uint32_t control_qubit_B_, uint32_t target_qubit_)
        : control_qubit_A{control_qubit_A_}, control_qubit_B{control_qubit_B_}, target_qubit{target_qubit_} {}
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t control_qubit_A;
    uint32_t control_qubit_B;
    uint32_t target_qubit;
};

class CircuitGate : public QGate
{
public:
    CircuitGate(const QCircuit& circuit_, uint32_t target_qubit_begin_, std::string name = "");
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override { return representation; }
    mutable QCircuit internal_circuit;
    std::string representation;
    uint32_t qubit_count;
    uint32_t target_qubit_begin;
};

class ControlCircuitGate : public QGate
{
public:
    ControlCircuitGate(const QCircuit& circuit_, uint32_t control_qubit_, uint32_t target_qubit_begin_, std::string name = "");
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override { return representation; }
    mutable QCircuit internal_circuit;
    std::string representation;
    uint32_t qubit_count;
    uint32_t control_qubit;
    uint32_t target_qubit_begin;
};


/**
 * @brief Prepares functionality for the ArrayFire Quantum Simulation functions
 * 
 * @details Initializes ArrayFire and arrays
 * 
 * @param argc 
 * @param argv 
 */
void initialize(int argc, char** argv);

//Single bit Operations
QState X_op(const QState& state);
QState Y_op(const QState& state);
QState Z_op(const QState& state);

QState RotateX_op(const QState& state, float angle);
QState RotateY_op(const QState& state, float angle);
QState RotateZ_op(const QState& state, float angle);

QState Hadamard_op(const QState& state);
QState Phase_op(const QState& state, float angle);

}