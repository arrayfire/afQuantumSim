/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h"

#include <random>
#include <stdexcept>
#include <unordered_map>

/*
    Pauli X matrix
    X = [ 0 , 1 ]
        [ 1 , 0 ]
*/
static const af::cfloat x_mat[4] = {
    {0.f, 0.f}, {1.f, 0.f}, {1.f, 0.f}, {0.f, 0.f}};

/*
    Pauli Y matrix
    Y = [ 0 , -i ]
        [ i ,  0 ]
*/
static const af::cfloat y_mat[4] = {
    {0.f, 0.f}, {0.f, -1.f}, {0.f, 1.f}, {0.f, 0.f}};

/*
    Pauli Z matrix
    Z = [ 1 ,  0 ]
        [ 0 , -1 ]
*/
static const af::cfloat z_mat[4] = {
    {1.f, 0.f}, {0.f, 0.f}, {0.f, 0.f}, {-1.f, 0.f}};

/*
    Hadamard matrix
    H = 1/sqrt(2) [ 1 ,  1 ]
                  [ 1 , -1 ]
*/
// clang format off
static const af::cfloat h_mat[4] = {
    {0.70710678118f , 0.f} , { 0.70710678118f , 0.f},
    {0.70710678118f , 0.f} , {-0.70710678118f , 0.f}
};
// clang format on

static af::array x_matrix;
static af::array y_matrix;
static af::array z_matrix;
static af::array hadamard_matrix;

static af::randomEngine intern_rnd_engine;

static std::random_device dv;
static std::mt19937 rd_gen{dv()};

static std::uniform_real_distribution<float> lin_dist{0.0f, 1.0f};

static std::unordered_map<std::string, std::shared_ptr<aqs::QCircuit>>
    cached_circuits;

namespace aqs {

void initialize(int argc, char** argv, af::Backend backend) {
    af::setBackend(backend);

    int device = (argc > 1) ? std::stoi(argv[1]) : 0;
    af::setDevice(device);

    x_matrix = af::array(2, 2, x_mat);

    y_matrix = af::array(2, 2, y_mat);

    z_matrix = af::array(2, 2, z_mat);

    hadamard_matrix = af::array(2, 2, h_mat);

    std::random_device rnd_device;
    intern_rnd_engine =
        af::randomEngine(AF_RANDOM_ENGINE_THREEFRY, rnd_device());
}

void clear_circuit_cache() { cached_circuits.clear(); }

QState::QState(const std::complex<float>& zeroState,
               const std::complex<float>& oneState)
    : state_{af::cfloat{zeroState.real(), zeroState.imag()},
             af::cfloat{oneState.real(), oneState.imag()}} {
    force_normalize();
}

QState::QState(const std::array<std::complex<float>, 2>& states)
    : state_{af::cfloat{states[0].real(), states[0].imag()},
             af::cfloat{states[1].real(), states[1].imag()}} {
    force_normalize();
}

QState& QState::set(const std::complex<float>& zero_state,
                    const std::complex<float>& one_state) {
    state_[0] = af::cfloat{zero_state.real(), zero_state.imag()};
    state_[1] = af::cfloat{one_state.real(), one_state.imag()};

    force_normalize();

    return *this;
}

bool QState::peek_measure() const {
    // Determine probability (value in range [0 , 1]) to obtain state |1>
    float prob = probability_true();

    // Generate randomly a value in range [0, 1]
    float val = lin_dist(rd_gen);

    // If random value is less than the |1> state probability return measurement
    // 1 Else return measurement 0
    return val < prob;
}

bool QState::measure() {
    // Get a random measurement
    bool measurement = peek_measure();

    // Update the state with the measurement
    state_[measurement]  = 1.f;
    state_[!measurement] = 0.f;

    return measurement;
}

std::array<uint32_t, 2> QState::profile_measure(uint32_t rep_count) const {
    uint32_t t = 0;
    for (uint32_t i = 0; i < rep_count; ++i)
        // Count the results that result in the |1> state
        t += static_cast<uint32_t>(peek_measure());

    return {rep_count - t, t};
}

void QState::force_normalize() {
    float mag2 =
        state_[0].real * state_[0].real + state_[0].imag * state_[0].imag +
        state_[1].real * state_[1].real + state_[1].imag * state_[1].imag;

    if (mag2 == 0.f)
        throw std::invalid_argument{"Cannot normalize a null state"};

    float mag = std::sqrt(mag2);

    state_[0] = state_[0] / mag;
    state_[1] = state_[1] / mag;
}

QCircuit::QCircuit(uint32_t qubit_count)
    : gate_list_{}
    , circuit_(
          af::identity(fast_pow2(qubit_count), fast_pow2(qubit_count), c32))
    , representation_{}
    , qubits_{qubit_count} {
    if (qubit_count < 1)
        throw std::invalid_argument{"Circuit must contain at least 1 qubit"};
    if (qubit_count > max_qubit_count)
        throw std::invalid_argument{"Maximum qubit count supported is " +
                                    std::to_string(max_qubit_count)};
}

void QCircuit::clear() {
    // Clear list of compiled gates
    clear_cache();

    // Reset list of added gates
    gate_list_.clear();
}

void QCircuit::clear_cache() {
    // Resets the cached index
    cached_index_ = 0;

    // Clear circuit matrix
    circuit_ = af::identity(state_count(), state_count(), c32);
}

bool operator==(const QCircuit& lhs, const QCircuit& rhs) {
    if (lhs.qubit_count() != rhs.qubit_count()) return false;
    if (lhs.gate_list().size() != rhs.gate_list().size()) return false;
    if (lhs.representation() != rhs.representation()) return false;
    for (uint32_t i = 0; i < lhs.gate_list().size(); ++i) {
        if (lhs.gate_list()[i] != rhs.gate_list()[i]) return false;
    }

    return true;
}

void QCircuit::compile() {
    if (cached_index_ != gate_list_.size()) {
        // Compiled non-cached gates
        for (std::size_t i = cached_index_; i < gate_list_.size(); ++i) {
            const auto& gate = *(gate_list_[i]);
            gate(*this);
        }

        // Update the last cached gate index
        cached_index_ = gate_list_.size();
    }
}

QSimulator::QSimulator(uint32_t qubit_count, const QState& initial_state,
                       const QNoise& noise_generator)
    : states_(qubit_count, initial_state)
    , statevector_(af::constant(0.f, fast_pow2(qubit_count), c32))
    , noise_{noise_generator}
    , qubits_{qubit_count}
    , basis_{Basis::Z} {
    if (initial_state == aqs::QState::zero())
        statevector_(0) = 1.f;
    else if (initial_state == aqs::QState::one())
        statevector_(state_count() - 1) = 1.f;
    else
        generate_statevector();
}

QSimulator::QSimulator(uint32_t qubit_count, std::vector<QState> initial_states,
                       const QNoise& noise_generator)
    : states_(std::move(initial_states))
    , statevector_(fast_pow2(qubit_count), c32)
    , noise_{noise_generator}
    , qubits_{qubit_count}
    , basis_{Basis::Z} {
    if (qubit_count != states_.size())
        throw std::invalid_argument{
            "The number of initial states must match the number of qubits in "
            "the circuit"};

    generate_statevector();
}

QSimulator::QSimulator(uint32_t qubit_count, const af::array& statevector,
                       const QNoise& noise_generator)
    : states_(qubit_count)
    , noise_{noise_generator}
    , qubits_{qubit_count}
    , basis_{Basis::Z} {
    if (statevector.dims()[0] != fast_pow2(qubit_count) ||
        statevector.dims()[1] != 1 || statevector.dims()[2] != 1 ||
        statevector.dims()[3] != 1)
        throw std::invalid_argument{"Invalid initial statevector shape"};

    auto norm = af::norm(statevector);

    if (norm == 0.)
        throw std::invalid_argument{"Cannot have a null statevector"};

    statevector_ = statevector / af::norm(statevector);
}

void QSimulator::generate_statevector() {
    std::vector<af::cfloat> temp(qubit_count() * 2);
    for (uint32_t i = 0; i < qubit_count(); ++i) {
        temp[i * 2]     = states_[i][0];
        temp[i * 2 + 1] = states_[i][1];
    }

    af::array states = af::array(2, qubit_count(), temp.data());

    statevector_ = states.col(0);

    // Generate the statevector from the tensor product of elements
    for (uint32_t i = 1; i < qubit_count(); ++i)
        statevector_ = tensor_product(statevector_, states.col(i));
}

void QSimulator::simulate(const QCircuit& circuit) {
    if (circuit.qubit_count() != qubits_)
        throw std::invalid_argument{
            "Number of qubit states and circuit input qubit states do not "
            "match"};

    statevector_ = af::matmul(circuit.circuit(), statevector_);

    QCircuit temp{qubit_count()};
    temp.circuit() = statevector_;
    for (std::size_t i = circuit.cached_index_; i < circuit.gate_list().size();
         ++i)
        circuit.gate_list()[i]->operator()(temp);
    statevector_ = temp.circuit();
}

bool QSimulator::peek_measure(uint32_t qubit) const {
    if (qubit >= qubit_count())
        throw std::out_of_range{"Cannot measure the state of the given qubit"};

    // Get probability of the qubit for state |1>
    float prob = qubit_probability_true(qubit);

    // Generate a random number
    float val = lin_dist(rd_gen);

    // If the random value is less than the probability |1> return 1
    // Else return 0
    return val < prob;
}

bool QSimulator::measure(uint32_t qubit) {
    if (qubit >= qubit_count())
        throw std::out_of_range{"Cannot measure the state of the given qubit"};

    // Generate a random number
    float val = lin_dist(rd_gen);

    // Find the probabilities of each state
    af::array probabilities = af::real(statevector_ * af::conjg(statevector_));

    // Mark the states at the position of the qubit
    af::array indices =
        af::iota(state_count(), 1, u32) & fast_pow2(qubit_count() - qubit - 1);

    // Select the states with qubit in |1> state
    af::array states1 = af::select(indices.as(b8), statevector_, 0.f);

    // Select the states with qubit in |0> state
    af::array states0 = af::select(indices.as(b8), 0.f, statevector_);

    // Find the probability for |1>
    float prob1 = af::sum<float>(states1 * af::conjg(states1));

    // If the random value is less than the probability |1> return 1
    // Else return 0
    bool measurement = val < prob1;

    // Update the statevector with the result of the qubit
    if (measurement)
        statevector_ = states1 / std::sqrt(prob1);
    else
        statevector_ = states0 / std::sqrt(1.f - prob1);

    return measurement;
}

uint32_t QSimulator::peek_measure_all() const {
    // Generate randomly a value in range [0, 1]
    float val = lin_dist(rd_gen);

    // Find the probabilities for each state
    af::array probabilities = af::real(statevector_ * af::conjg(statevector_));

    // Find all the states that posses a probability greater than or equal to
    // the random value generated
    af::array prob = af::where(af::ceil(af::accum(probabilities) - val));

    // The entry of that first value is the state measured
    uint32_t state = prob.isempty() ? 0 : prob(0).scalar<uint32_t>();

    return state;
}

uint32_t QSimulator::measure_all() {
    // Generate random measurement
    uint32_t measurement = peek_measure_all();

    // Set the statevector to the result of the measurement
    statevector_(af::span)    = 0.f;
    statevector_(measurement) = 1.f;

    return measurement;
}

float QSimulator::qubit_probability_true(uint32_t qubit) const {
    if (qubit >= qubit_count() || qubit < 0)
        throw std::out_of_range{"Cannot obtain probability of the given qubit"};

    // Find the probabilities of each state
    af::array probabilities = af::real(statevector_ * af::conjg(statevector_));

    // Mark the states at the position of the qubit
    af::array indices =
        af::iota(state_count(), 1, u32) & fast_pow2(qubit_count() - qubit - 1);

    // Select the states with qubit in |1> state
    af::array states1 = af::select(indices.as(b8), statevector_, 0.);

    // Sum the probabilities of all the states that contain the qubit in state
    // |1>
    float prob1 = af::sum<float>(states1 * af::conjg(states1));

    return prob1;
}

float QSimulator::state_probability(uint32_t state) const {
    if (state >= state_count() || state < 0)
        throw std::out_of_range{"Cannot obtain probability of the given state"};

    // Find the state
    af::cfloat val = statevector_(state).scalar<af::cfloat>();

    // Return probability (state magnitude)
    return val.real * val.real + val.imag * val.imag;
}

std::vector<float> QSimulator::probabilities() const {
    std::vector<float> out(state_count());

    // Find the probabilities of each state
    af::array probs = af::real(statevector_ * af::conjg(statevector_));

    // Copy to vector
    probs.host(out.data());

    return out;
}

void QSimulator::set_basis(Basis basis) {
    if (basis != basis_) {
        // Hadamard X
        // clang format off
        static af::cfloat z_to_x_vals[] = {
            { 0.70710678118f, 0.f } , { 0.70710678118f , 0.f },
            { 0.70710678118f, 0.f } , {-0.70710678118f , 0.f }
        };
        // clang format on

        // Hadamard Y
        // clang format off
        static af::cfloat z_to_y_vals[] = {
            { 0.70710678118f, 0.f } , { 0.70710678118f , 0.f },
            { 0.f, -0.70710678118f } , { 0.f , 0.70710678118f }
        };
        // clang format on

        static af::array z_to_x = af::array(2, 2, z_to_x_vals).T();
        static af::array z_to_y = af::array(2, 2, z_to_y_vals).T();
        static af::array x_to_z = af::inverse(z_to_x);
        static af::array y_to_z = af::inverse(z_to_y);

        // Tranform from current basis to z-basis
        af::array matrix;
        switch (basis_) {
            case Basis::Z: matrix = af::identity(2, 2, c32); break;
            case Basis::Y: matrix = y_to_z;
            case Basis::X: matrix = x_to_z; break;
        }

        auto gen_tensor_matrix = [qubits = qubits_](const af::array& matrix) {
            af::array out = matrix;
            for (uint32_t i = 0; i < qubits - 1; ++i)
                out = tensor_product(matrix, out);
            return out;
        };

        // Transform from z-basis to target basis
        switch (basis) {
            case Basis::Z: break;
            case Basis::Y: matrix = af::matmul(z_to_y, matrix); break;
            case Basis::X: matrix = af::matmul(z_to_x, matrix); break;
        }

        // Update statevector
        statevector_ = af::matmul(gen_tensor_matrix(matrix), statevector_);
        basis_       = basis;
    }
}

std::vector<uint32_t> QSimulator::profile_measure_all(
    uint32_t rep_count) const {
    const auto states = state_count();
    std::vector<uint32_t> count(states);

    // Generate rep_count amount of random numbers in [0, 1] and repeat them
    // along dim 1
    af::array rnd =
        af::tile(af::randu(rep_count, f32, intern_rnd_engine), 1, states);

    // Generate the cumulative probability of each state along dim 1 and repeat
    // them along dim 0
    af::array probabilities = af::tile(
        af::transpose(
            af::accum(af::real(statevector_ * af::conjg(statevector_)))),
        rep_count, 1);

    // Find all the states that posses a probability greater than or equal to
    // the random value generated as indices + 1
    auto temp = af::sum(af::ceil(probabilities - rnd), 1).as(u32);

    // Sort the indices of the states measured
    af::array prob =
        af::sort(af::join(0, temp, af::iota(states, 1, u32) + 1), 0, false);

    // Count the indices of the states measured
    af::array keys, vals;
    af::countByKey(keys, vals, prob, af::constant(1, rep_count + states));
    vals -= 1;

    // Copy it to the host
    vals.host(count.data());

    return count;
}

std::array<uint32_t, 2> QSimulator::profile_measure(uint32_t qubit,
                                                    uint32_t rep_count) const {
    if (qubit >= qubit_count())
        throw std::out_of_range{
            "Cannot profile measurement of the given qubit"};

    auto states         = state_count();
    int32_t qubit_state = 1 << (qubit_count() - qubit - 1);

    // Generate rep_count of random floats between 0 and 1
    af::array rnd = af::randu(rep_count, f32, intern_rnd_engine);

    // Find the probability of all states
    af::array probabilities = af::real(statevector_ * af::conjg(statevector_));

    // Mark the states where the qubit is in state 1
    af::array one_indices = (af::iota(states, 1, u32) & qubit_state).as(b8);

    // Obtain the probability of the states with qubit in state 1
    af::array one_vals = af::select(one_indices, probabilities, 0LL);

    // Obtain the probability of being in state 1 and tile it with rep_count
    af::array one_sum = af::tile(af::sum(one_vals), rep_count);

    // Find all the random numbers which result in state |1> and count them
    auto temp = af::sum<uint32_t>(af::ceil(one_sum - rnd), 0);

    return {rep_count - temp, temp};
}

bool X::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& X::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);

    // Flip the bits of the position of the target qubit
    auto col_indices = af::iota(states, 1, s32) ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    // Generate pauli-X sparse matrix
    auto matrix_x =
        af::sparse(states, states, af::constant(af::cfloat{1.f, 0.f}, states),
                   row_indices, col_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_x, circuit);

    return qc;
}

bool X::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const X*>(&rhs)->target_qubit;
}

std::string X::to_string() const {
    return "X,0,1:" + std::to_string(target_qubit) + ";";
}

bool Y::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& Y::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota          = af::iota(states, 1, s32);

    // Flip the bits of the position of the target qubit
    auto col_indices = iota ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    // Set all |0> positions to i
    auto values = af::constant(af::cfloat{0.f, 1.f}, states);

    // Replace all |1> positions to -i
    af::replace(values, (iota & (1 << (qubits - target_qubit - 1))).as(b8),
                af::constant(af::cfloat{0.f, -1.f}, states));

    // Generate pauli-Y sparse matrix
    auto matrix_y =
        af::sparse(states, states, values, row_indices, col_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_y, circuit);

    return qc;
}

bool Y::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const Y*>(&rhs)->target_qubit;
}

std::string Y::to_string() const {
    return "Y,0,1:" + std::to_string(target_qubit) + ";";
}

bool Z::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& Z::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota          = af::iota(states, 1, s32);
    auto col_indices   = iota;
    auto row_indices   = af::iota(states + 1, 1, s32);

    // Set all |1> positions to -1
    auto vals = af::constant(af::cfloat{-1.f, 0.f}, states);

    // Replace all |0> positions to 1
    af::replace(vals, (iota & mask).as(b8),
                af::constant(af::cfloat{1.f, 0.f}, states));

    // Generate pauli-Z sparse matrix
    auto matrix_z = af::sparse(states, states, vals, row_indices, col_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_z, circuit);

    return qc;
}

bool Z::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const Z*>(&rhs)->target_qubit;
}

std::string Z::to_string() const {
    return "Z,0,1:" + std::to_string(target_qubit) + ";";
}

bool RotX::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& RotX::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity =
        af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity =
        af::identity(fast_pow2(qubits - target_qubit - 1),
                     fast_pow2(qubits - target_qubit - 1), c32);

    // Create rotation X matrix
    auto cos_angle          = std::cos(angle / 2.0f);
    auto sin_angle          = std::sin(angle / 2.0f);
    const af::cfloat vals[] = {{cos_angle, 0.f},
                               {0.f, -sin_angle},
                               {0.f, -sin_angle},
                               {cos_angle, 0.f}};

    // Find the n-qubit RotX matrix for the target qubit as I_L @ RotX @ I_R
    // where @ is the tensor product
    af::array temp = tensor_product(left_identity, af::array(2, 2, vals).T());
    temp           = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit       = af::matmul(temp, circuit);

    return qc;
}

bool RotX::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const RotX*>(&rhs)->target_qubit &&
           angle == static_cast<const RotX*>(&rhs)->angle;
}

std::string RotX::to_string() const {
    return "RotX,0,1:" + std::to_string(target_qubit) + ";";
}

bool RotY::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& RotY::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity =
        af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity =
        af::identity(fast_pow2(qubits - target_qubit - 1),
                     fast_pow2(qubits - target_qubit - 1), c32);

    // Create rotation Y matrix
    auto cos_angle          = std::cos(angle / 2.0f);
    auto sin_angle          = std::sin(angle / 2.0f);
    const af::cfloat vals[] = {{cos_angle, 0.f},
                               {-sin_angle, 0.f},
                               {sin_angle, 0.f},
                               {cos_angle, 0.f}};

    // Find the n-qubit RotY matrix for the target qubit as I_L @ RotY @ I_R
    // where @ is the tensor product
    af::array temp = tensor_product(left_identity, af::array(2, 2, vals).T());
    temp           = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit       = af::matmul(temp, circuit);

    return qc;
}

bool RotY::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const RotY*>(&rhs)->target_qubit &&
           angle == static_cast<const RotY*>(&rhs)->angle;
}

std::string RotY::to_string() const {
    return "RotY,0,1:" + std::to_string(target_qubit) + ";";
}

bool RotZ::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& RotZ::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity =
        af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity =
        af::identity(fast_pow2(qubits - target_qubit - 1),
                     fast_pow2(qubits - target_qubit - 1), c32);

    // Create rotation Z matrix
    auto cos_angle          = std::cos(angle / 2.0f);
    auto sin_angle          = std::sin(angle / 2.0f);
    const af::cfloat vals[] = {{cos_angle, -sin_angle},
                               {0.f, 0.f},
                               {0.f, 0.f},
                               {cos_angle, sin_angle}};

    // Find the n-qubit RotZ matrix for the target qubit as I_L @ RotZ @ I_R
    // where @ is the tensor product
    af::array temp = tensor_product(left_identity, af::array(2, 2, vals));
    temp           = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit       = af::matmul(temp, circuit);

    return qc;
}

bool RotZ::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const RotZ*>(&rhs)->target_qubit &&
           angle == static_cast<const RotZ*>(&rhs)->angle;
}

std::string RotZ::to_string() const {
    return "RotZ,0,1:" + std::to_string(target_qubit) + ";";
}

bool H::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& H::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity =
        af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity =
        af::identity(fast_pow2(qubits - target_qubit - 1),
                     fast_pow2(qubits - target_qubit - 1), c32);

    // Find the n-qubit X matrix for the target qubit as I_L @ H @ I_R where @
    // is the tensor product
    af::array temp = tensor_product(left_identity, hadamard_matrix.T());
    temp           = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit       = af::matmul(temp, circuit);
    return qc;
}

bool H::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const H*>(&rhs)->target_qubit;
}

std::string H::to_string() const {
    return "H,0,1:" + std::to_string(target_qubit) + ";";
}

bool Phase::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    return true;
}

QCircuit& Phase::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota          = af::iota(states, 1, s32);
    auto col_indices   = iota;
    auto row_indices   = af::iota(states + 1, 1, s32);

    // Generate the values for the phase part
    auto vals =
        af::constant(af::cfloat{std::cos(angle), std::sin(angle)}, states);

    // Replace |1> with the phase rotation
    af::replace(vals, (iota & mask).as(b8),
                af::constant(af::cfloat{1.f, 0.f}, states));

    // Generate phase matrix
    auto matrix_phase =
        af::sparse(states, states, vals, row_indices, col_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_phase, circuit);

    return qc;
}

bool Phase::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const Phase*>(&rhs)->target_qubit &&
           angle == static_cast<const Phase*>(&rhs)->angle;
}

std::string Phase::to_string() const {
    if (angle == pi / 2)
        return "S,0,1:" + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 2)
        return "S†,0,1:" + std::to_string(target_qubit) + ";";
    else if (angle == pi / 4)
        return "T,0,1:" + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 4)
        return "T†,0,1:" + std::to_string(target_qubit) + ";";
    else
        return "Phase,0,1:" + std::to_string(target_qubit) + ";";
}

bool Swap::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{
            "Cannot use the swap gate on the same target qubits"};
    return true;
}

QCircuit& Swap::operator()(QCircuit& qc) const {
    const auto states = qc.state_count();
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{
            "Cannot use the swap gate on the same target qubits"};

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int32_t maskA = 1 << (qubits - 1 - target_qubit_A);
    int32_t maskB = 1 << (qubits - 1 - target_qubit_B);

    // Create a mask for flipping the bits
    auto col_indices = af::iota(states, 1, s32);
    auto temp        = ((col_indices >> posA) ^ (col_indices >> posB)) & 1;

    // Generate column indices
    // Swap the positions of the column by switching the bits in the qubit
    // position
    col_indices = col_indices ^ ((temp << posA) | (temp << posB));

    auto row_indices = af::iota(states + 1, 1, s32);

    // Generate swap sparse matrix
    auto swap_matrix =
        af::sparse(states, states, af::constant(af::cfloat{1.0f, 0.f}, states),
                   row_indices, col_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(swap_matrix, circuit);

    return qc;
}

bool Swap::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit_A == static_cast<const Swap*>(&rhs)->target_qubit_A &&
           target_qubit_B == static_cast<const Swap*>(&rhs)->target_qubit_B;
}

std::string Swap::to_string() const {
    return "Swap,0,2:" + std::to_string(target_qubit_A) + "," +
           std::to_string(target_qubit_B) + ";";
}

bool CX::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    return true;
}

QCircuit& CX::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask  = 1 << (qubits - 1 - target_qubit);

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Flip the bits of the position of the target qubit if the control qubit is
    // set to |1>
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8),
                iota ^ target_mask);

    // Generate CX sparse matrix
    auto cx_matrix =
        af::sparse(states, states, af::constant(af::cfloat{1.f, 0.f}, states),
                   row_indices, column_indices);

    // Update the circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(cx_matrix, circuit);

    return qc;
}

bool CX::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CX*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CX*>(&rhs)->control_qubit;
}

std::string CX::to_string() const {
    return "X,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CY::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    return true;
}

QCircuit& CY::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask  = 1 << (qubits - 1 - target_qubit);
    int32_t mask         = control_mask | target_mask;

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Flip the bits of the position of the target qubit if the control qubit is
    // set to |1>
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8),
                iota ^ target_mask);

    // Set all states to be the same by default
    auto values     = af::constant(af::cfloat{1.f, 0.f}, states);
    auto op_indices = iota & mask;

    // Replace the |0> with i if the control qubit is set to |1>
    af::replace(values, op_indices != mask,
                af::constant(af::cfloat{0.f, 1.f}, states));

    // Replace the |1> with -i if the control qubit is set to |1>
    af::replace(values, op_indices != control_mask,
                af::constant(af::cfloat{0.f, -1.f}, states));

    // Generate CY sparse matrix
    auto cy_matrix =
        af::sparse(states, states, values, row_indices, column_indices);

    // Update the circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(cy_matrix, circuit);

    return qc;
}

bool CY::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CY*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CY*>(&rhs)->control_qubit;
}

std::string CY::to_string() const {
    return "Y,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CZ::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    return true;
}

QCircuit& CZ::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask  = 1 << (qubits - 1 - target_qubit);
    int32_t mask         = control_mask | target_mask;

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Set all states be set to 1.0 by default
    auto values = af::constant(af::cfloat{1.f, 0.f}, states);

    // Replace all |1> states if the control qubit is set to |1>
    af::replace(values, (iota & mask) != mask,
                af::constant(af::cfloat{-1.f, 0.f}, states));

    // Generate CZ sparse matrix
    auto cz_matrix =
        af::sparse(states, states, values, row_indices, column_indices);

    // Update the circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(cz_matrix, circuit);

    return qc;
}

bool CZ::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CZ*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CZ*>(&rhs)->control_qubit;
}

std::string CZ::to_string() const {
    return "Z,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CPhase::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    return true;
}

QCircuit& CPhase::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask  = 1 << (qubits - 1 - target_qubit);
    int32_t mask         = control_mask | target_mask;

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Set all states by default to 1.0f
    auto values = af::constant(af::cfloat{1.f, 0.f}, states);

    // Replace all states containing |1> with the phase rotation if the control
    // qubit is set to |1>
    af::replace(
        values, (iota & mask) != mask,
        af::constant(af::cfloat{std::cos(angle), std::sin(angle)}, states));

    // Generate CPhase sparse matrix
    auto cphase_matrix =
        af::sparse(states, states, values, row_indices, column_indices);

    // Update the circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(cphase_matrix, circuit);

    return qc;
}

bool CPhase::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CPhase*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CPhase*>(&rhs)->control_qubit &&
           angle == static_cast<const CPhase*>(&rhs)->angle;
}

std::string CPhase::to_string() const {
    if (angle == pi / 2)
        return "S,1,1:" + std::to_string(control_qubit) + "," +
               std::to_string(target_qubit) + ";";
    else if (angle == -pi / 2)
        return "S†,1,1:" + std::to_string(control_qubit) + "," +
               std::to_string(target_qubit) + ";";
    else if (angle == pi / 4)
        return "T,1,1:" + std::to_string(control_qubit) + "," +
               std::to_string(target_qubit) + ";";
    else if (angle == -pi / 4)
        return "T†,1,1:" + std::to_string(control_qubit) + "," +
               std::to_string(target_qubit) + ";";
    else
        return "Phase,1,1:" + std::to_string(control_qubit) + "," +
               std::to_string(target_qubit) + ";";
}

bool CSwap::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given circuit"};
    if (target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit_A || control_qubit == target_qubit_B)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{
            "Cannot use the swap gate on the same target qubits"};
    return true;
}

QCircuit& CSwap::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given circuit"};
    if (target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit_A || control_qubit == target_qubit_B)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{
            "Cannot use the swap gate on the same target qubits"};

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int32_t control = 1 << (qubits - control_qubit - 1);
    int32_t targetA = 1 << posA;
    int32_t targetB = 1 << posB;

    auto iota        = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);

    // Generate column indices that should be swapped (swap the bits in the
    // states of qA and qB)
    auto temp         = (iota >> posA) ^ (iota >> posB);
    temp              = temp & 1;
    auto replace_vals = iota ^ ((temp << posA) | (temp << posB));

    auto column_indices = iota;
    af::replace(column_indices, (iota & control) != control, replace_vals);

    // Generate CSwap sparse matrix
    auto matrix_cswap =
        af::sparse(states, states, af::constant(af::cfloat{1.f, 0.f}, states),
                   row_indices, column_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_cswap, circuit);

    return qc;
}

bool CSwap::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit_A == static_cast<const CSwap*>(&rhs)->target_qubit_A &&
           target_qubit_B == static_cast<const CSwap*>(&rhs)->target_qubit_B &&
           control_qubit == static_cast<const CSwap*>(&rhs)->control_qubit;
}

std::string CSwap::to_string() const {
    return "Swap,1,2:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit_A) + std::to_string(target_qubit_B) +
           ";";
}

bool CH::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};
    return true;
}

QCircuit& CH::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    const int32_t value_mask         = (1 << (qubits - target_qubit)) | 1;
    const int32_t control_index_mask = 1 << (qubits - control_qubit);
    const int32_t index_mask         = ~(1 << (qubits - target_qubit - 1));

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32) * 2;
    auto resized_iota   = af::flat(af::tile(iota.T(), 2));
    auto column_indices = resized_iota & index_mask;

    const int32_t off[] = {0, 1 << (qubits - target_qubit - 1)};

    auto offset = af::tile(af::array(2, off), states);
    column_indices += offset;

    auto iota2               = af::iota(states * 2, 1, s32);
    auto control_indices     = (iota2 & control_index_mask).as(b8);
    auto new_col_indices     = af::shift(resized_iota, -1);
    new_col_indices(af::end) = states;

    af::replace(column_indices, control_indices, new_col_indices);

    const float sqrt2      = 0.70710678118f;
    const int32_t val_mask = (1 << (qubits - target_qubit)) | 1;

    const af::cfloat identity_vals[] = {{1.f, 0.f}, {0.f, 0.f}};

    // Select all positions where the control qubit is |1>
    auto values        = af::tile(af::array(2, identity_vals), states);
    auto value_indices = ((iota2 & value_mask) ^ value_mask).as(b8);

    // Replace the enabled places with the hadamard matrix
    af::replace(values, !(control_indices && value_indices),
                af::constant(af::cfloat{sqrt2, 0.f}, states * 2));
    af::replace(values, !(control_indices && !value_indices),
                af::constant(af::cfloat{-sqrt2, 0.f}, states * 2));

    // Generate sparse matrix
    auto matrix_chadamard =
        af::sparse(states, states, values, row_indices, column_indices);

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(matrix_chadamard, circuit);

    return qc;
}

bool CH::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CH*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CH*>(&rhs)->control_qubit;
}

std::string CH::to_string() const {
    return "H,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CRotX::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    return true;
}

QCircuit& CRotX::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    // Apply a Control Gate with RotX as the target gate
    ControlGate(RotX::gate(angle), control_qubit, target_qubit)(qc);

    return qc;
}

bool CRotX::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CRotX*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CRotX*>(&rhs)->control_qubit &&
           angle == static_cast<const CRotX*>(&rhs)->angle;
}

std::string CRotX::to_string() const {
    return "RotX,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CRotY::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    return true;
}

QCircuit& CRotY::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    // Apply a Control Gate with RotY as the target gate
    ControlGate(RotY::gate(angle), control_qubit, target_qubit)(qc);

    return qc;
}

bool CRotY::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CRotY*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CRotY*>(&rhs)->control_qubit &&
           angle == static_cast<const CRotY*>(&rhs)->angle;
}

std::string CRotY::to_string() const {
    return "RotY,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CRotZ::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    return true;
}

QCircuit& CRotZ::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit == control_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    // Apply a Control Gate with RotZ as the target gate
    ControlGate(RotZ::gate(angle), control_qubit, target_qubit)(qc);

    return qc;
}

bool CRotZ::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CRotZ*>(&rhs)->target_qubit &&
           control_qubit == static_cast<const CRotZ*>(&rhs)->control_qubit &&
           angle == static_cast<const CRotZ*>(&rhs)->angle;
}

std::string CRotZ::to_string() const {
    return "RotZ,1,1:" + std::to_string(control_qubit) + "," +
           std::to_string(target_qubit) + ";";
}

bool CCNot::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();

    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_A == target_qubit || control_qubit_B == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    return true;
}

QCircuit& CCNot::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_A == target_qubit || control_qubit_B == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = (1 << (qubits - 1 - control_qubit_A)) |
                           (1 << (qubits - 1 - control_qubit_B));
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Replace all places where both control qubits are |1> with swapped indices
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8),
                iota ^ target_mask);

    // Create CCNot sparse matrix
    auto ccnot_matrix =
        af::sparse(states, states, af::constant(af::cfloat{1.f, 0.f}, states),
                   row_indices, column_indices);

    // Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit       = af::matmul(ccnot_matrix, circuit);

    return qc;
}

bool CCNot::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const CCNot*>(&rhs)->target_qubit &&
           control_qubit_A ==
               static_cast<const CCNot*>(&rhs)->control_qubit_A &&
           control_qubit_B == static_cast<const CCNot*>(&rhs)->control_qubit_B;
}

std::string CCNot::to_string() const {
    return "X,2,1:" + std::to_string(control_qubit_A) + "," +
           std::to_string(control_qubit_B) + "," +
           std::to_string(target_qubit) + ";";
}

bool Or::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();

    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_A == target_qubit || control_qubit_B == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    return true;
}

QCircuit& Or::operator()(QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_A == target_qubit || control_qubit_B == target_qubit)
        throw std::invalid_argument{
            "Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = (1 << (qubits - 1 - control_qubit_A)) |
                           (1 << (qubits - 1 - control_qubit_B));
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);

    auto iota           = af::iota(states, 1, s32);
    auto row_indices    = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    // Replace all indices where control qubits or result in |1> with swapped
    // indices
    af::replace(column_indices, ((iota & control_mask) == 0).as(b8),
                iota ^ target_mask);

    // Generate Or sparse matrix
    auto or_matrix =
        af::sparse(states, states, af::constant(af::cfloat{1.f, 0.f}, states),
                   row_indices, column_indices);

    // Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit       = af::matmul(or_matrix, circuit);

    return qc;
}

bool Or::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit == static_cast<const Or*>(&rhs)->target_qubit &&
           control_qubit_A == static_cast<const Or*>(&rhs)->control_qubit_A &&
           control_qubit_B == static_cast<const Or*>(&rhs)->control_qubit_B;
}

std::string Or::to_string() const {
    std::stringstream buffer;

    // Equivalent circuit for 2-qubit Or with 1 target
    buffer << "P;"
           << "X,0,1:" << control_qubit_A << ";"
           << "X,0,1:" << control_qubit_B << ";"
           << "X,0,1:" << target_qubit << ";"
           << "X,2,1:" << control_qubit_A << "," << control_qubit_B << ","
           << target_qubit << ";"
           << "X,0,1:" << control_qubit_A << ";"
           << "X,0,1:" << control_qubit_B << ";"
           << "P;";

    return buffer.str();
}

/**
 * @brief Given the circuit string representation and starting target qubit,
 *        it updates the string to reflect the new positions for the internal
 * gates added
 *
 * @param circuit_string circuit string representation to be updated
 * @param offset index of the target qubit begin
 *
 * @return std::string circuit string representation for the updated circuit
 */
static std::string update_circuit_representation(
    const std::string& circuit_string, uint32_t offset);

/**
 * @brief Given the circuit string representation, starting target qubit, and
 * additional control qubit; it updates the string to reflect the new positions
 * for the internal gates added
 *
 * @param circuit_string circuit string representation to be updated
 * @param control index of the control qubit
 * @param target index of the target qubit begin
 * @return std::string circuit string representation for the updated circuit
 */
static std::string update_ctrl_circuit_representation(
    const std::string& circuit_string, uint32_t control, uint32_t target);

Gate::Gate(const QCircuit& circuit_, uint32_t target_qubit_begin_,
           std::string name)
    : representation{}
    , qubit_count{circuit_.qubit_count()}
    , target_qubit_begin{target_qubit_begin_} {
    auto circuit_iter = cached_circuits.find(circuit_.representation());
    if (circuit_iter == cached_circuits.end()) {
        auto res = cached_circuits.insert(
            {circuit_.representation(), std::make_shared<QCircuit>(circuit_)});
        if (!res.second)
            throw std::runtime_error{"Could not store circuit"};
        else
            circuit_iter = res.first;
    } else {
        if (!(circuit_ == *(circuit_iter->second)))
            cached_circuits[circuit_.representation()] =
                std::make_shared<QCircuit>(circuit_);
        circuit_iter = cached_circuits.find(circuit_.representation());
    }
    internal_circuit = circuit_iter->second;

    if (name == "") {
        representation = update_circuit_representation(
            circuit_.representation(), target_qubit_begin_);
    } else {
        std::stringstream buff;

        // Naming of the circuit must not contain parsing tokens
        if (name.find_first_of(",;:") != std::string::npos)
            throw std::invalid_argument{
                "Name cannot contain commas, colons, nor semicolons"};

        // Set circuit to have no ctrl qubits
        buff << name << ",0," << circuit_.qubit_count() << ":";

        // Offset all the target qubits by target_qubit_begin
        buff << target_qubit_begin_;
        for (uint32_t i = 1; i < circuit_.qubit_count(); ++i)
            buff << "," << i + target_qubit_begin_;
        buff << ";";

        representation = buff.str();
    }
}

bool Gate::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit_begin >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    return true;
}

QCircuit& Gate::operator()(QCircuit& qc) const {
    const auto qubits      = qc.qubit_count();
    const auto states      = qc.state_count();
    const auto state_count = fast_pow2(qubit_count);
    if (target_qubit_begin >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    auto& circuit = qc.circuit();
    internal_circuit->compile();

    if (qc.qubit_count() != internal_circuit->qubit_count()) {
        const uint32_t circuit_qubits   = qc.qubit_count();
        const uint32_t circuit_states   = qc.state_count();
        const uint32_t gate_qubits      = qubit_count;
        const uint32_t gate_states      = 1 << gate_qubits;
        const uint32_t gate_qubit_begin = target_qubit_begin;
        af::array gate_matrix =
            af::identity(circuit_states, circuit_states, c32);

        // Generate the matrix of the source circuit
        const af::array& gate = internal_circuit->circuit();

        uint32_t rem_count = 1 << (circuit_qubits - gate_qubits);

        auto len = gate_states * gate_states * rem_count;
        auto m   = af::tile(
              af::flat(af::tile(af::iota(gate_states, 1, s32).T(), gate_states)),
              rem_count);
        auto n   = af::iota(gate_states, gate_states * rem_count, s32);
        auto ind = af::flat(af::tile(af::iota(rem_count, 1, s32).T(),
                                     gate_states * gate_states));

        // Find the target indices where the values of the source circuit should
        // go
        auto ii = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, n,
                            ind, len);
        auto jj = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, m,
                            ind, len);

        // Get a list of the entries in order where should they be placed
        af::array gate_values = gate(n * gate_states + m);

        // Replace the entries of the target matrix with the entries selected
        // from the source matrix
        gate_matrix(ii * circuit_states + jj) = gate_values;

        // Update circuit
        circuit = af::matmul(gate_matrix, circuit);
    } else
        circuit = af::matmul(internal_circuit->circuit(), circuit);

    return qc;
}

bool Gate::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit_begin ==
               static_cast<const Gate*>(&rhs)->target_qubit_begin &&
           internal_circuit == static_cast<const Gate*>(&rhs)->internal_circuit;
}

ControlGate::ControlGate(const QCircuit& circuit_, uint32_t control_qubit_,
                         uint32_t target_qubit_begin_, std::string name)
    : representation{}
    , qubit_count{circuit_.qubit_count()}
    , control_qubit{control_qubit_}
    , target_qubit_begin{target_qubit_begin_} {
    auto circuit_iter = cached_circuits.find(circuit_.representation());
    if (circuit_iter == cached_circuits.end()) {
        auto res = cached_circuits.insert(
            {circuit_.representation(), std::make_shared<QCircuit>(circuit_)});
        if (!res.second)
            throw std::runtime_error{"Could not store circuit"};
        else
            circuit_iter = res.first;
    } else {
        if (!(circuit_ == *(circuit_iter->second)))
            cached_circuits[circuit_.representation()] =
                std::make_shared<QCircuit>(circuit_);
        circuit_iter = cached_circuits.find(circuit_.representation());
    }
    internal_circuit = circuit_iter->second;

    if (name == "") {
        representation = update_ctrl_circuit_representation(
            circuit_.representation(), control_qubit, target_qubit_begin_);
    } else {
        std::stringstream buff;

        // Naming of the circuit must not contain parsing tokens
        if (name.find_first_of(",;:") != std::string::npos)
            throw std::invalid_argument{
                "Name cannot contain commas, colons, nor semicolons"};

        // Set the circuit to contain 1 ctrl qubit and append it to the front of
        // the qubit list
        buff << name << ",1," << circuit_.qubit_count() << ":" << control_qubit_
             << ",";

        // Offset all the target qubits by target_qubit_begin
        buff << target_qubit_begin_;
        for (uint32_t i = 1; i < circuit_.qubit_count(); ++i)
            buff << "," << i + target_qubit_begin_;
        buff << ";";

        representation = buff.str();
    }
}

bool ControlGate::check(const QCircuit& qc) const {
    const auto qubits = qc.qubit_count();
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Gate must fit inside the circuit qubit count"};
    if (target_qubit_begin <= control_qubit &&
        control_qubit < target_qubit_begin + qubit_count)
        throw std::out_of_range{
            "Control qubit cannot be one of the target qubits of the gate"};
    if (qubit_count >= qubits)
        throw std::invalid_argument{"Cannot add a bigger gate to the circuit"};
    if (control_qubit >= qubits)
        throw std::out_of_range{
            "Control qubit must be inside the circuit qubit range"};

    return true;
}

QCircuit& ControlGate::operator()(QCircuit& qc) const {
    const auto qubits      = qc.qubit_count();
    const auto states      = qc.state_count();
    const auto state_count = fast_pow2(qubit_count);

    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Gate must fit inside the circuit qubit count"};
    if (target_qubit_begin <= control_qubit &&
        control_qubit < target_qubit_begin + qubit_count)
        throw std::out_of_range{
            "Control qubit cannot be one of the target qubits of the gate"};
    if (qubit_count + 1 > qubits)
        throw std::invalid_argument{"Cannot add a bigger gate to the circuit"};

    const uint32_t circuit_qubits   = qc.qubit_count();
    const uint32_t circuit_states   = qc.state_count();
    const uint32_t gate_qubits      = qubit_count;
    const uint32_t gate_states      = 1 << gate_qubits;
    const uint32_t gate_qubit_begin = target_qubit_begin;
    af::array gate_matrix = af::identity(circuit_states, circuit_states, c32);

    // Generate the matrix of the source circuit
    internal_circuit->compile();
    const af::array& gate = internal_circuit->circuit();

    uint32_t rem_count = 1 << (circuit_qubits - gate_qubits - 1);
    auto len           = gate_states * gate_states * rem_count;
    auto m             = af::tile(
                    af::flat(af::tile(af::iota(gate_states, 1, s32).T(), gate_states)),
                    rem_count);
    auto n   = af::iota(gate_states, gate_states * rem_count, s32);
    auto ind = af::flat(
        af::tile(af::iota(rem_count, 1, s32).T(), gate_states * gate_states));
    auto offset = control_qubit < target_qubit_begin ? 0 : 1;

    // Find the target indices where the values of the source circuit should go
    auto ii = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits - offset,
                        n, ind, len);
    auto jj = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits - offset,
                        m, ind, len);

    // The target indices must be those where the control qubit is |1>
    auto ones = af::constant(1, len, s32);
    ii        = insert_bits(ii, ones,
                            af::constant(qubits - control_qubit - 1, len, s32), ones,
                            len);
    jj        = insert_bits(jj, ones,
                            af::constant(qubits - control_qubit - 1, len, s32), ones,
                            len);

    // Get a list of the entries in order where should they be placed
    af::array gate_values = gate(n * gate_states + m);

    // Replace the entries of the target matrix with the entries selected from
    // the source matrix
    gate_matrix(ii * circuit_states + jj) = gate_values;

    // Update circuit
    auto& circuit = qc.circuit();
    circuit       = af::matmul(gate_matrix, circuit);

    return qc;
}

bool ControlGate::operator==(const QGate& rhs) const noexcept {
    return type() == rhs.type() &&
           target_qubit_begin ==
               static_cast<const ControlGate*>(&rhs)->target_qubit_begin &&
           control_qubit ==
               static_cast<const ControlGate*>(&rhs)->control_qubit &&
           internal_circuit ==
               static_cast<const ControlGate*>(&rhs)->internal_circuit;
}

QState X_op(const QState& state) {
    af::cfloat temp[2];
    temp[0] = (x_mat[0] * state[0]) + (x_mat[1] * state[1]);
    temp[1] = (x_mat[2] * state[0]) + (x_mat[3] * state[1]);

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState Y_op(const QState& state) {
    af::cfloat temp[2];
    temp[0] = (y_mat[0] * state[0]) + (y_mat[1] * state[1]);
    temp[1] = (y_mat[2] * state[0]) + (y_mat[3] * state[1]);

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState Z_op(const QState& state) {
    af::cfloat temp[2];
    temp[0] = (z_mat[0] * state[0]) + (z_mat[1] * state[1]);
    temp[1] = (z_mat[2] * state[0]) + (z_mat[3] * state[1]);

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState RotateX_op(const QState& state, float angle) {
    af::cfloat temp[2];
    temp[0] = state[0] * std::cos(angle / 2.f) +
              state[1] * af::cfloat{0.f, -std::sin(angle / 2.f)};
    temp[1] = state[0] * af::cfloat{0.f, -std::sin(angle / 2.f)} +
              state[1] * std::cos(angle / 2.f);

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState RotateY_op(const QState& state, float angle) {
    af::cfloat temp[2];
    temp[0] =
        state[0] * std::cos(angle / 2.f) + state[1] * -std::sin(angle / 2.f);
    temp[1] =
        state[0] * std::sin(angle / 2.f) + state[1] * std::cos(angle / 2.f);

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState RotateZ_op(const QState& state, float angle) {
    af::cfloat temp[2];
    temp[0] =
        state[0] * af::cfloat{std::cos(angle / 2.f), -std::sin(angle / 2.f)} +
        state[1] * 0.f;
    temp[1] = state[0] * 0.f + state[1] * af::cfloat{std::cos(angle / 2.f),
                                                     std::sin(angle / 2.f)};

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState Hadamard_op(const QState& state) {
    af::cfloat temp[2];
    temp[0] = h_mat[0] * state[0] + h_mat[1] * state[1];
    temp[1] = h_mat[2] * state[0] + h_mat[3] * state[1];

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

QState Phase_op(const QState& state, float angle) {
    af::cfloat temp[2];
    temp[0] = state[0] * 1.f + state[1] * 0.f;
    temp[1] = state[0] * 0.f +
              state[1] * af::cfloat{std::cos(angle), std::sin(angle)};

    return {{temp[0].real, temp[0].imag}, {temp[1].real, temp[1].imag}};
}

std::string update_circuit_representation(const std::string& circuit_string,
                                          uint32_t offset) {
    std::stringstream out;
    std::size_t current_begin = 0;
    auto current_end          = circuit_string.find(";", current_begin);

    while (current_end != std::string::npos) {
        auto name_end = circuit_string.find(",", current_begin);
        auto name_str =
            circuit_string.substr(current_begin, name_end - current_begin);
        if (name_end > current_end) {
            out << circuit_string.substr(current_begin,
                                         current_end - current_begin)
                << ";";
            current_begin = current_end + 1;
            current_end   = circuit_string.find(";", current_begin);
            continue;
        }

        out << name_str << ",";

        auto colon_pos       = circuit_string.find(":", current_begin);
        auto qubit_pos_begin = name_end + 1;
        if (colon_pos < current_end) {
            out << circuit_string.substr(name_end + 1, colon_pos - name_end);
            qubit_pos_begin = colon_pos + 1;
        }

        auto qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        qubit_pos_end =
            qubit_pos_end > current_end ? current_end : qubit_pos_end;

        // Add the offset to all target qubits
        while (qubit_pos_end < current_end) {
            out << std::stoi(circuit_string.substr(
                       qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) +
                       offset
                << ",";
            qubit_pos_begin = qubit_pos_end + 1;
            qubit_pos_end   = circuit_string.find(",", qubit_pos_begin);
        }
        qubit_pos_end =
            qubit_pos_end > current_end ? current_end : qubit_pos_end;
        out << std::stoi(circuit_string.substr(
                   qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) +
                   offset
            << ";";

        current_begin = current_end + 1;
        current_end   = circuit_string.find(";", current_begin);
    }

    return out.str();
}

std::string update_ctrl_circuit_representation(
    const std::string& circuit_string, uint32_t control, uint32_t target) {
    std::stringstream out;
    std::size_t current_begin = 0;
    auto current_end          = circuit_string.find(";", current_begin);

    while (current_end != std::string::npos) {
        auto name_end = circuit_string.find(",", current_begin);
        auto name_str =
            circuit_string.substr(current_begin, name_end - current_begin);
        if (name_end > current_end) {
            out << circuit_string.substr(current_begin,
                                         current_end - current_begin)
                << ";";
            current_begin = current_end + 1;
            current_end   = circuit_string.find(";", current_begin);
            continue;
        }
        out << name_str << ",";

        // Add ctrl to all qubits
        auto colon_pos       = circuit_string.find(":", current_begin);
        auto qubit_pos_begin = name_end + 1;
        if (colon_pos < current_end) {
            auto temp = circuit_string.find(",", name_end + 1);
            out << std::stoi(circuit_string.substr(name_end + 1,
                                                   temp - name_end - 1)) +
                       1
                << circuit_string.substr(temp, colon_pos - temp + 1);
            qubit_pos_begin = colon_pos + 1;
        }

        out << control << ",";

        auto qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        qubit_pos_end =
            qubit_pos_end > current_end ? current_end : qubit_pos_end;

        // Add the offset to all target qubits
        while (qubit_pos_end < current_end) {
            out << std::stoi(circuit_string.substr(
                       qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) +
                       target
                << ",";
            qubit_pos_begin = qubit_pos_end + 1;
            qubit_pos_end   = circuit_string.find(",", qubit_pos_begin);
        }
        qubit_pos_end =
            qubit_pos_end > current_end ? current_end : qubit_pos_end;
        out << std::stoi(circuit_string.substr(
                   qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) +
                   target
            << ";";

        current_begin = current_end + 1;
        current_end   = circuit_string.find(";", current_begin);
    }

    return out.str();
}

}  // namespace aqs