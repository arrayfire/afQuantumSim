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
 * 
 * ******* Arrayfire Quantum Simulator AQS **********
 * 
 * ////// How to use //////
 * 1. Initialize the library by calling `aqs::initialize()` and passing the command line arguments to it
 * 
 *    e.g. int main(int argc, char** argv) {
 *              ...
 *              aqs::initialize(argc, argv);
 *              ...
 *         }
 * 
 * 
 * 2. Create a Quantum Circuit by using the class `aqs::QCircuit`.
 *    Specify the number of qubits by passing it to the contructor
 * 
 *    e.g. aqs::QCircuit qc{2};
 * 
 * 3. Insert gates by using the `<<` (left shift) operator with the Quantum Circuit at the left and any added gates at the right.
 *    To create the gates, call the constructor of each gate. Depending on the gate being added,
 *    you may need to specify to which qubit the gate wil be added.
 * 
 *    e.g. qc << aqs::H{0} << aqs::CNOT{0, 1};
 * 
 * 4. Call the QCircuit member function `aqs::QCircuit::compile` to generate the internal matrix representation
 *    of the circuit.
 *    Note that this function can be call anywhere, but it will only generate the matrix of the all gates
 *    that have been added prior to that point. This function executes the compilation (math calculations) 
 *    of the circuit, so it is performance intensive; for more than 10 qubits it might take a while
 * 
 *    e.g. qc.compile();
 * 
 * 5. To simulate, create the context of the quantum simulator with the class `aqs::QSimulator`
 *    It must be initialized with the same number of qubits as the circuit,
 *    and can be added an initial state to setup the statevector in.
 * 
 *    e.g. aqs::QSimulator qs{2, aqs::QState::one()}; // All the qubits are set into an inital |1> state
 * 
 * 6. Finally, to execute the simulation, call the QSimulator Member function `aqs::QSimulator::simulate`
 *    using the simulator created and pass the circuit as the argument.
 *    This will result in a complete Quantum Computer Simulation
 * 
 *    Note that this also a performance intensive function, and the result of the simulation is stored in the QSimultor
 * 
 *    e.g. qs.simulate(qc);
 * 
 * 
 * 7. To obtain information from the simulation, you may use `aqs::QSimulator::statevector` to obtain the statevector,
 *    call `aqs::QSimulator::measure` to collapse the state of a given qubit and return the output;
 *    `aqs::QSimulator::measure_all` to collapse the state of all the qubits and return it as binary string;
 *    `aqs::QSimulator::peek_measure` and `aqs::QSimulator::peek_measure_all()` for the non-collapsing/non-affecting versions for it;
 *    and `aqs::QSimulator::profile_measure` and `aqs::QSimulator::profile_measure_all` to execute a random measurements multiple times
 * 
 *    e.g. auto result = qs.peek_measure_all(); // get one random result for all the qubit measured states
 *         auto profile = qs.profile_measure(2, 100); // measure 100 times the collpased state of the simulation of qubit 2
 * 
 * 
 * ////// Important Implementation Details //////
 *  - The internal simulation of a Quantum Computer's state is done using a statevector instead of a density matrix
 *  - The qubits are indexed from 0 to qubit_count - 1.
 *  - The first qubit (qubit 0) is the most significant qubit and the last qubit is least signficant qubit in terms of positioning
 *    This means that first qubit is always first in the tensor product and the last qubit will go last in the tensor proudct,
 *    and that measuring will result in the result of the last qubit being stored in the least significant bit of the binary string
 *  - The representation of complex states and calculations use single floating point numbers.
 *  - The maximum number of qubits that can be used in the simulations is 30.
 *  - All the matrix and vectors are stored in ArrayFire arrays of type c32 (1 complex float = 2 floats)
 * 
 */


/**
 * @brief ArrayFire Quantum Simulator
 * 
 */
namespace aqs
{

static constexpr float pi = 3.14159265358979323846f;
static constexpr uint32_t max_qubit_count = 30;

class QGate;
class QState;
class QCircuit;
class QSimulator;

/**
 * @brief Prepares functionality for the ArrayFire Quantum Simulation functions
 * 
 * @details Initializes ArrayFire and internal arrays
 * 
 * @param argc 
 * @param argv 
 * @param backend arrayfire backend for the simulator to use
 */
void initialize(int argc, char** argv, af::Backend backend = af::Backend::AF_BACKEND_DEFAULT);

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


/**
 * @brief Class managing the behavior of a Quantum Circuit
 *        including the addition of gates, the circuit matrix,
 *        and its representation
 * 
 */
class QCircuit
{
public:
    friend class QSimulator;

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

    template<typename T>
    friend QCircuit& operator<<(QCircuit& qc, const T& gate);

    template<typename T>
    friend QCircuit& operator<<(QCircuit& qc, const std::vector<T>& gates);

    /**
     * @brief Returns the number of qubit registers in the circuit
     * 
     * @return uint32_t 
     */
    uint32_t qubit_count() const noexcept { return qubits_; }

    /**
     * @brief Returns the number of different unique orthogonal states
     *        that the state vector can be represented in
     * 
     * @return uint32_t 
     */
    uint32_t state_count() const noexcept { return fast_pow2(qubits_); }

    /**
     * @brief Returns a reference to the current matrix representation of the circuit
     * 
     * @warning If the array is modified, it must maintain its unitary property and dimensions 
     *          Modifying the internal circuit incorrectly may result in undefined behavior
     *          Only use for low level access and functionality that the gates/circuit do not provide
     * 
     * @return af::array& 
     */
    af::array& circuit() noexcept { return circuit_; }

    /**
     * @brief Returns a const reference to the current matrix representation of the circuit
     * 
     * @return const af::array& 
     */
    const af::array& circuit() const noexcept { return circuit_; }

    /**
     * @brief Returns a reference to the list of gates that have been added
     * 
     * @return std::vector<std::shared_ptr<aqs::QGate>>&
     */
    auto& gate_list() noexcept { return gate_list_; }

    /**
     * @brief Returns a const reference to the list of gates that have been added
     * 
     * @return const std::vector<std::shared_ptr<aqs::QGate>>&
     */
    const auto& gate_list() const noexcept { return gate_list_; }

    std::string& representation() noexcept { return representation_; }

    const std::string& representation() const noexcept { return representation_; }

    /**
     * @brief Starts the computation of the matrix representation of the circuit
     *        from the added gates
     * 
     * @note Adding gates after this call will not update matrix representation,
     *       a new call for this function must be added to update it
     * 
     * @details This function caches its result to compute only the gates that 
     *          have not been computed into the circuit matrix
     * 
     */
    void compile();

    /**
     * @brief Clears the circuit by removing all gates and resetting the internal circuit matrix
     * 
     */
    void clear();

private:
    std::vector<std::shared_ptr<QGate>> gate_list_;
    af::array circuit_;
    std::string representation_;
    uint32_t qubits_ = 0;
    std::size_t cached_index_ = 0;
};

/**
 * @brief Class managing the logic for the simulation of quantum noise models
 * 
 */
class QNoise
{
};

/**
 * @brief  Class managing the simulation of a quantum circuit and a initial state
 * 
 */
class QSimulator
{
public:

    /**
     * @brief Enum of the different types Basis that the qubits can be measured in
     * 
     */
    enum class Basis : int8_t
    {
        // |0_z> = |0>
        // |1_z> = |1>
        Z,

        // |0_y> = 1/sqrt(2)|0> + i/sqrt(2) |1>
        // |1_y> = 1/sqrt(2)|0> - i/sqrt(2) |1>
        Y,

        // |0_x> = 1/sqrt(2)|0> + 1/sqrt(2) |1>
        // |1_x> = 1/sqrt(2)|0> - 1/sqrt(2) |1>
        X
    };

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
     * @brief Construct a new QSimulator object
     * 
     * @param qubit_count number of qubits of the simulator
     * @param statevector the statevector initialized to
     * @param noise_generator the noise generator to be used
     */
    explicit QSimulator(uint32_t qubit_count,
                        const af::array& statevector,
                        const QNoise& noise_generator = QNoise{});

    /**
     * @brief Calculates the statevector from all the individual states of the qubits
     * 
     */
    void generate_statevector();
    
    /**
     * @brief Computes the simulation for the circuit with the given initial conditions
     *        and updates the statevector to the result
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
     * @brief Measures the given qubit and collapses its state (and updates the statevector) to it
     * 
     * @details In the statevector, all states that contain the given measurement are updated according to the previous probability
     *          using conditional probability, while the other states are set to 0
     * 
     * @param qubit index of the qubit
     * 
     * @return true for |1> state measurement
     * @return false for |0> state measurement
     */
    bool measure(uint32_t qubit);

    /**
     * @brief Returns a bit representation of the statevector measured after measuring all bits.
     *        It collapses the statevector to that state
     * 
     * @return uint32_t measurement
     */
    uint32_t measure_all();

    /**
     * @brief Returns a bit representation of the statevector measured after measuring all bits.
     *        It does not collpase the statevector
     * 
     * @return uint32_t measurement
     */
    uint32_t peek_measure_all() const;

    /**
     * @brief Profiles the measurement of the given state of a qubit from the stored statevector using the given number of tests
     * 
     * @param qubit qubit to measure
     * @param rep_count number of measurements to be done for profiling
     * @return std::array<int, 2> [0] = # of measurements for state |0>; [1] = # of measurements for state |1>
     */
    std::array<uint32_t, 2> profile_measure(uint32_t qubit, uint32_t rep_count) const;

    /**
     * @brief Profiles the measurement of the given output state from the stored statevector using the given number of tests
     * 
     * @param rep_count number of measurements to be done for profiling
     * @return std::vector<int> vector with the measurements [int(xxxx...)] = #number of measurements of state |xxxxx> where xxxxx is the binary output state
     */
    std::vector<uint32_t> profile_measure_all(uint32_t rep_count) const;

    /**
     * @brief Returns a reference to the initial QState of given qubit from the simulator
     * 
     * @param index index of the qubit
     * @return QState& 
     */
    QState& qubit(uint32_t index) noexcept { assert(index < qubit_count()); return states_[index]; };

    /**
     * @brief Returns a reference to the initial QState of given qubit from the simulator
     * 
     * @param index index of the qubit
     * @return const QState& 
     */
    const QState& qubit(uint32_t index) const noexcept { assert(index < qubit_count()); return states_[index]; };

    /**
     * @brief Returns the complex number of the state in the statevector vector
     * 
     * @param state int representation of the binary unique state
     * @return af::cfloat 
     */
    af::cfloat state(uint32_t state) const noexcept { assert(state < state_count()); return statevector_(state).scalar<af::cfloat>(); }

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
     * @brief Returns the probability to measure the given state from teh statevector
     * 
     * @param state int representation of the binary unique state
     * @return float probability of |xxxxxx> where xxxxx is the state
     */
    float state_probability(uint32_t state) const;

    /**
     * @brief Returns a list of the probabilities for all the states
     * 
     * @return std::vector<float> 
     */
    std::vector<float> probabilities() const;

    /**
     * @brief Returns the number of qubits of the simulator
     * 
     * @return int number of qubits
     */
    uint32_t qubit_count() const noexcept { return qubits_; }

    /**
     * @brief Returns the number of unique states that can be measured from the statevector
     * 
     * @return int number of states
     */
    uint32_t state_count() const noexcept { return fast_pow2(qubits_);}

    /**
     * @brief Returns the basis meaasurement direction
     * 
     * @return Basis 
     */
    Basis get_basis() const noexcept { return basis_; }

    /**
     * @brief Changes the basis for measurement to the passed direction
     * 
     * @param basis basis direction
     */
    void set_basis(Basis basis);

    /**
     * @brief Returns a af::array& to the internal statevector stored by the simulator
     * 
     * @warning If the array is modified, it must maintain its normalized property and dimensions 
     *          Modifying the internal statevector incorrectly may result in undefined behavior
     *          Only use for low level access and functionality that the simulator do not provide
     * 
     * @return af::array& 
     */
    af::array& statevector() noexcept { return statevector_; }

    /**
     * @brief Returns a const af::array& to the internal statevector stored by the simulator
     * 
     * @return const af::array& 
     */
    const af::array& statevector() const noexcept { return statevector_; }

private:
    std::vector<QState> states_;
    af::array statevector_;
    QNoise noise_;
    uint32_t qubits_;
    Basis basis_ = Basis::Z;
};


/**
 * @brief Pure virtual base class delineating the structure for Quantum Gates
 * 
 */
class QGate
{
public:
    /**
     * @brief Applies the gate matrix to the circuit matrix
     * 
     * @param qc circuit to apply the gate to
     * @return QCircuit& reference to the input circuit
     */
    virtual QCircuit& operator()(QCircuit& qc) const = 0;

    /**
     * @brief Returns a string representation of the gate
     * 
     * @return std::string 
     */
    virtual std::string to_string() const = 0;

    /**
     * @brief Returns a unsigned integer that uniquely identifies
     *        the gate
     * 
     * @return uint32_t 
     */
    virtual uint32_t type() const noexcept = 0;

    /**
     * @brief Checks if the gate can be added to the passed circuit
     * 
     * @note This function may throw with the error before returning false
     * 
     * @param qc circuit in which the gate will be added
     * @return true the gate will be added successfully to the circuit
     * @return false the gate will not be added to the circuit
     */
    virtual bool check(const QCircuit& qc) const = 0;

protected:

    /**
     * @brief Lists of unique identifiers for the internal gates implemented
     * 
     */
    enum class GateTypes : uint32_t
    {
        Barrier,
        X, Y, Z, Hadamard, Phase, Swap, RotX, RotY, RotZ,
        CX, CY, CZ, CHadamard, CPhase, CSwap, CRotX, CRotY, CRotZ,
        CCX, Or, Circuit, ControlCircuit
    };
};

template<typename T>
QCircuit& operator<<(QCircuit& qc, const T& gate)
{
    static_assert(std::is_base_of<QGate, T>::value, "Gate must inherit from QGate class");
    if (gate.check(qc))
    {
        qc.representation_.append(gate.to_string());
        qc.gate_list().push_back(std::make_shared<T>(gate));
    }

    return qc;
}

template<typename T>
QCircuit& operator<<(QCircuit& qc, const std::vector<T>& gates)
{
    static_assert(std::is_base_of<QGate, T>::value, "Gate must inherit from QGate class");
    for (const auto& gate : gates)
    {
        if (gate.check(qc))
        {
            qc.representation_.append(gate.to_string());
            qc.gate_list().push_back(std::make_shared<T>(gate));
        }
    }

    return qc;
}

/**
 * @brief Barrier gate: marks a separation between prior and subsequent gates in the circuit
 *        Can be made visible or invisible when displaying the circuit
 * 
 */
class Barrier : public QGate
{
public:
    Barrier(bool visible_ = true) noexcept : visible{visible_} {}

    bool check(const QCircuit&) const override { return true; }
    QCircuit& operator()(QCircuit& qc) const override { return qc; }
    std::string to_string() const override { return visible ? "B;" : "P;"; }
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Barrier); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Barrier); }

    bool visible = true;
};


/**
 * @brief Pauli X gate: applies the pauli X matrix on the given qubit
 *        (Equivalent to applying a NOT gate)
 * 
 */
class X : public QGate
{
public:
    X(uint32_t target_qubit_) noexcept : target_qubit{target_qubit_} {}

    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::X); }
    bool check(const QCircuit&) const override;

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::X); }
    static const QCircuit& gate() {
        static QCircuit qc = [](){ QCircuit qc(1); qc << X{0}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
};

/**
 * @brief NOT gate: applies a NOT gate to the given qubit
 *        (Equivalent to applying X gate)
 * 
 */
using Not = X;

/**
 * @brief Pauli Y gate: applies the pauli Y gate matrix on the given qubit
 * 
 */
class Y : public QGate
{
public:
    Y(uint32_t target_qubit_) noexcept : target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Y); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Y); }
    static const QCircuit& gate() {
        static QCircuit qc = [](){ QCircuit qc(1); qc << Y{0}; qc.compile(); return qc; }();
        return qc;
    }
 
    uint32_t target_qubit;
};

/**
 * @brief Pauli Z gate: applies the Pauli Z matrix to the given qubit
 *        (Equivalent to a Phase-Pi gate)
 * 
 */
class Z : public QGate
{
public:
    Z(uint32_t target_qubit_) noexcept : target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Z); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Z); }
    static const QCircuit& gate() {
        static QCircuit qc = [](){ QCircuit qc(1); qc << Z{0}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
};

/**
 * @brief Rotation X gate: applies a rotation around the x-axis to the given qubit 
 * 
 */
class RotX : public QGate
{
public:
    RotX(uint32_t target_qubit_, float angle_) noexcept : target_qubit{target_qubit_} , angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::RotX); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::RotX); }
    static QCircuit gate(float angle) {
        QCircuit qc = [&angle](){ QCircuit qc(1); qc << RotX{0, angle}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Rotation Y gate: applies a rotation around the y-axis to the given qubit 
 * 
 */
class RotY : public QGate
{
public:
    RotY(uint32_t target_qubit_, float angle_) noexcept : target_qubit{target_qubit_} , angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::RotY); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::RotY); }
    static QCircuit gate(float angle) {
        QCircuit qc = [&angle](){ QCircuit qc(1); qc << RotY{0, angle}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Rotation Z gate: applies a rotation around the z-axis to the given qubit 
 * 
 */
class RotZ : public QGate
{
public:
    RotZ(uint32_t target_qubit_, float angle_) noexcept : target_qubit{target_qubit_} , angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::RotZ); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::RotZ); }
    static QCircuit gate(float angle) {
        QCircuit qc = [&angle](){ QCircuit qc(1); qc << RotZ{0, angle}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Hadamard gate: applies a Hadamard gate (superposition gate) to the given qubit
 *        (Equivalent to transforming the state from the Z-basis states to X-basis states)
 * 
 */
class H : public QGate
{
public:
    H(uint32_t target_qubit_) noexcept : target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Hadamard); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Hadamard); }
    static const QCircuit& gate() {
        static QCircuit qc = [](){ QCircuit qc(1); qc << H{0}; qc.compile(); return qc; }();
        return qc;
    }

    uint32_t target_qubit;
};

/**
 * @brief Phase gate: Applies a phase rotation to the given qubit
 *        (Equivalent to rotating the |1> component around in the xy-plane)
 * 
 */
class Phase : public QGate
{
public:
    Phase(uint32_t target_qubit_, float angle_) noexcept : target_qubit{target_qubit_}, angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Phase); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Phase); }
    static QCircuit gate(float angle){
        return [angle](){ QCircuit qc{ 1 }; qc << Phase{ 0, angle }; qc.compile(); return qc; }();
    }

    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Swap gate: Swaps the state of two qubits
 * 
 */
class Swap : public QGate
{
public:
    Swap(uint32_t target_qubit_A_, uint32_t target_qubit_B_) noexcept
        : target_qubit_A{target_qubit_A_} , target_qubit_B{target_qubit_B_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Swap); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Swap); }

    uint32_t target_qubit_A;
    uint32_t target_qubit_B;
};

/**
 * @brief Control Pauli X gate: Applies Pauli X gate on a target qubit
 *        depending on the state of a control qubit
 *        (Equivalent to a CNOT/XOR gate)
 */
class CX : public QGate
{
public:
    CX(uint32_t control_qubit_, uint32_t target_qubit_) noexcept
        : control_qubit{control_qubit_} , target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CX); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CX); }

    uint32_t control_qubit;
    uint32_t target_qubit;
};

/**
 * @brief Control Pauli Y gate: Applies Pauli Y gate on a target qubit
 *        depending on the state of a control qubit
 * 
 */
class CY : public QGate
{
public:
    CY(uint32_t control_qubit_, uint32_t target_qubit_) noexcept
        : control_qubit{control_qubit_} , target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CY); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CY); }

    uint32_t control_qubit;
    uint32_t target_qubit;
};

/**
 * @brief Control Pauli Z gate: Applies Pauli Z gate on a target qubit
 *        depending on the state of a control qubit
 * 
 */
class CZ : public QGate
{
public:
    CZ(uint32_t control_qubit_, uint32_t target_qubit_) noexcept
        : control_qubit{control_qubit_} , target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CZ); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CZ); }

    uint32_t control_qubit;
    uint32_t target_qubit;
};

/**
 * @brief Control Hadamard gate: Applies Hadamard gate on a target qubit
 *        depending on the state of a control qubit
 * 
 */
class CH : public QGate
{
public:
    CH(uint32_t control_qubit_, uint32_t target_qubit_) noexcept
        : control_qubit{control_qubit_} , target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CHadamard); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CHadamard); }

    uint32_t control_qubit;
    uint32_t target_qubit;
};

/**
 * @brief Control Not gate: Applies NOT gate on a target qubit
 *        depending on the state of a control qubit
 *        (Equivalent to a CX/XOR gate)
 * 
 */
using CNot = CX;

/**
 * @brief Logical XOR gate: Returns the result XOR of two qubits in the target qubit
 * 
 */
using Xor = CX;

/**
 * @brief Control Phase gate: Applies a phase rotation to the target qubit
 *        depending on the state of a control qubit
 * 
 */
class CPhase : public QGate
{
public:
    CPhase(uint32_t control_qubit_, uint32_t target_qubit_, float angle_)
        : control_qubit{control_qubit_} , target_qubit{target_qubit_} , angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CPhase); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CPhase); }

    uint32_t control_qubit;
    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Control Swap gate: Swaps the state of two qubits depending on the state of a control qubit
 * 
 */
class CSwap : public QGate
{
public:
    CSwap(uint32_t control_qubit_, uint32_t target_qubit_A_, uint32_t target_qubit_B_) noexcept
        : control_qubit{control_qubit_}, target_qubit_A{target_qubit_A_}, target_qubit_B{target_qubit_B_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CSwap); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CSwap); }

    uint32_t control_qubit;
    uint32_t target_qubit_A;
    uint32_t target_qubit_B;
};

/**
 * @brief Control Rotation X gate: applies a rotation around the x-axis to the given qubit
 *        depending on the state of a control qubit
 * 
 */
class CRotX : public QGate
{
public:
    CRotX(uint32_t control_qubit_, uint32_t target_qubit_, float angle_) noexcept
        : control_qubit{ control_qubit_ } , target_qubit{ target_qubit_ } , angle{angle_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CRotX); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CRotX); }

    uint32_t control_qubit;
    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Control Rotation Y gate: applies a rotation around the y-axis to the given qubit
 *        depending on the state of a control qubit
 * 
 */
class CRotY : public QGate
{
public:
    CRotY(uint32_t control_qubit_, uint32_t target_qubit_, float angle_) noexcept
        : control_qubit{ control_qubit_ } , target_qubit{ target_qubit_ } , angle{ angle_ } {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CRotY); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CRotY); }

    uint32_t control_qubit;
    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Control Rotation Z gate: applies a rotation around the z-axis to the given qubit
 *        depending on the state of a control qubit
 * 
 */
class CRotZ : public QGate
{
public:
    CRotZ(uint32_t control_qubit_, uint32_t target_qubit_, float angle_) noexcept
        : control_qubit{ control_qubit_ } , target_qubit{ target_qubit_ } , angle{ angle_ } {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CRotZ); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CRotZ); }

    uint32_t control_qubit;
    uint32_t target_qubit;
    float angle;
};

/**
 * @brief Control-Control Not gate: A double-qubit controlled not gate
 *        Applies a NOT gate depending on the state of two control qubits
 *        (Equivalent to a CCX)
 * 
 */
class CCNot : public QGate
{
public:
    CCNot(uint32_t control_qubit_A_, uint32_t control_qubit_B_, uint32_t target_qubit_)
        : control_qubit_A{control_qubit_A_}, control_qubit_B{control_qubit_B_}, target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::CCX); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::CCX); }

    uint32_t control_qubit_A;
    uint32_t control_qubit_B;
    uint32_t target_qubit;
};

/**
 * @brief Logical AND gate: Returns the result AND of control two qubits in the target qubit
 * 
 */
using And = CCNot;

/**
 * @brief Logical OR gate: Returns the result OR of control two qubits in the target qubit
 * 
 */
class Or : public QGate
{
public:
    Or(uint32_t control_qubit_A_, uint32_t control_qubit_B_, uint32_t target_qubit_) noexcept
        : control_qubit_A{control_qubit_A_}, control_qubit_B{control_qubit_B_}, target_qubit{target_qubit_} {}

    std::string to_string() const override;
    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Or); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Or); }

    uint32_t control_qubit_A;
    uint32_t control_qubit_B;
    uint32_t target_qubit;
};

/**
 * @brief Circuit Gate: Allows for the insertion of Quantum Circuits as a custom gate
 * 
 */
class Gate : public QGate
{
public:

    /**
     * @brief Constructs a Circuit Gate
     * 
     * @param circuit_ the circuit to insert as a gate
     * @param target_qubit_begin_ beginning of the range
     * @param name alias for the gate (empty string leaves the internal circuit visible)
     */
    Gate(const QCircuit& circuit_, uint32_t target_qubit_begin_, std::string name = "");

    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override { return representation; }
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::Circuit); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::Circuit); }

    mutable QCircuit internal_circuit;
    std::string representation;
    uint32_t qubit_count;
    uint32_t target_qubit_begin;
};


/**
 * @brief Control Circuit gate: Allows for the insertion of a controllable Quantum Circuit as a gate
 * 
 */
class ControlGate : public QGate
{
public:
    /**
     * @brief Construct a Control Circuit Gate
     * 
     * @param circuit_ circuit to insert as a gate
     * @param control_qubit_ control qubit for the gate
     * @param target_qubit_begin_ beginning of the range of the gate 
     * @param name name alias for the gate (empty string leaves the internal circuit visible)
     */
    ControlGate(const QCircuit& circuit_, uint32_t control_qubit_, uint32_t target_qubit_begin_, std::string name = "");

    bool check(const QCircuit&) const override;
    QCircuit& operator()(QCircuit&) const override;
    std::string to_string() const override { return representation; }
    uint32_t type() const noexcept override { return static_cast<uint32_t>(GateTypes::ControlCircuit); }

    static constexpr uint32_t static_type() noexcept { return static_cast<uint32_t>(GateTypes::ControlCircuit); }

    mutable QCircuit internal_circuit;
    std::string representation;
    uint32_t qubit_count;
    uint32_t control_qubit;
    uint32_t target_qubit_begin;
};

/**
 * @brief Returns pauli X matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState X_op(const QState& state);

/**
 * @brief Returns pauli Y matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState Y_op(const QState& state);

/**
 * @brief Returns pauli Z matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState Z_op(const QState& state);

/**
 * @brief Returns Rotate X matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState RotateX_op(const QState& state, float angle);

/**
 * @brief Returns Rotate Y matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState RotateY_op(const QState& state, float angle);

/**
 * @brief Returns Rotate Z matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState RotateZ_op(const QState& state, float angle);

/**
 * @brief Returns Hadamard matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState Hadamard_op(const QState& state);

/**
 * @brief Returns Phase matrix applied on the given qubit
 * 
 * @param state state to be applied on
 * @return QState 
 */
QState Phase_op(const QState& state, float angle);

}