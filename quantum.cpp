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
    {0.f , 0.f} , {1.f , 0.f},
    {1.f , 0.f} , {0.f , 0.f}
};

/*
    Pauli Y matrix
    Y = [ 0 , -i ]
        [ i ,  0 ] 
*/
static const af::cfloat y_mat[4] = {
    {0.f , 0.f} , {0.f , -1.f},
    {0.f , 1.f} , {0.f ,  0.f}
};

/*
    Pauli Z matrix
    Z = [ 1 ,  0 ]
        [ 0 , -1 ]
*/
static const af::cfloat z_mat[4] = {
    {1.f , 0.f} , {  0.f , 0.f},
    {0.f , 0.f} , { -1.f , 0.f}
};

/*
    Hadamard matrix
    H = 1/sqrt(2) [ 1 ,  1 ]
                  [ 1 , -1 ]
*/
static const af::cfloat h_mat[4] = {
    {0.70710678118f , 0.f} , { 0.70710678118f , 0.f},
    {0.70710678118f , 0.f} , {-0.70710678118f , 0.f}
};

static af::array x_matrix;
static af::array y_matrix;
static af::array z_matrix;
static af::array hadamard_matrix;

static af::randomEngine intern_rnd_engine;

static std::random_device dv;
static std::mt19937 rd_gen{dv()};

static std::uniform_real_distribution<float> lin_dist(0.0f, 1.0f);
static std::unordered_map<uint64_t, af::array> cached_arrays;

namespace aqs
{

enum class cache_array_type : int
{
    iota, 
};

static uint64_t gen_cache_array_index(cache_array_type arr_type, int32_t size)
{
    return (static_cast<int64_t>(arr_type) << 32) | size;
}

static af::array gen_array(uint64_t identifier)
{
    int32_t size = static_cast<int32_t>(identifier & static_cast<uint64_t>(~0));
    cache_array_type type = static_cast<cache_array_type>(identifier >> 32);

    switch (type)
    {
    case cache_array_type::iota:
        return af::iota(size, 1, s32);
    default:
        break;
    }
}

static const af::array& get_cached_array(uint64_t identifier)
{
    auto iter = cached_arrays.find(identifier);
    if (iter == cached_arrays.end())
        return cached_arrays[identifier] = gen_array(identifier);
    return iter->second;
}

QState::QState(const std::complex<float>& zeroState, const std::complex<float>& oneState)
    :state_{ af::cfloat{zeroState.real(), zeroState.imag()} , af::cfloat{oneState.real(), oneState.imag()}}
{
    force_normalize();
}

QState::QState(const std::array<std::complex<float>, 2>& states)
    :state_{af::cfloat{states[0].real() , states[0].imag()}, af::cfloat{states[1].real() , states[1].imag()}}
{   
    force_normalize();
}

QState& QState::set(const std::complex<float>& zero_state, const std::complex<float>& one_state)
{
    state_[0] = af::cfloat{zero_state.real() , zero_state.imag()};
    state_[1] = af::cfloat{ one_state.real() , one_state.imag()};

    force_normalize();

    return *this;
}

bool QState::peek_measure() const
{
    // Determine probability (value in range [0 , 1]) to obtain state |1>
    float prob = probability_true();

    // Generate randomly a value in range [0, 1]
    float val = lin_dist(rd_gen);

    // If random value is greater than the |1> state probability return measurement 1
    // Else return measurement 0
    return val < prob;
}

bool QState::measure()
{
    bool measurement = peek_measure();

    state_[ measurement] = 1.f;
    state_[!measurement] = 0.f;

    return measurement;
}

std::array<uint32_t, 2> QState::profile_measure(uint32_t rep_count) const
{
    uint32_t t = 0;
    for (uint32_t i = 0; i < rep_count; ++i)
        t += static_cast<uint32_t>(peek_measure());

    return { rep_count - t, t };
}

void QState::force_normalize()
{
    float mag2 = state_[0].real * state_[0].real + state_[0].imag * state_[0].imag +
                 state_[1].real * state_[1].real + state_[1].imag * state_[1].imag;
    
    if (mag2 == 0.f)
        throw std::invalid_argument{"Cannot normalize a null state"};
    
    float mag = std::sqrt(mag2);
    
    state_[0] = state_[0] / mag;
    state_[1] = state_[1] / mag;
}

QCircuit::QCircuit(uint32_t qubit_count)
    :gate_list_{}, circuit_(af::identity(fast_pow2(qubit_count), fast_pow2(qubit_count), c32)), representation_{}, qubits_{qubit_count}
{
    if (qubit_count < 1)
        throw std::invalid_argument{"Circuit must contain at least 1 qubit"};
    if (qubit_count > max_qubit_count)
        throw std::invalid_argument{"Maximum qubit count supported is " + std::to_string(max_qubit_count)};
}
    
void QCircuit::Global_Measure()
{
}

void QCircuit::Measure(uint32_t qubit)
{
    const auto count = qubit_count();
    if (qubit >= count || qubit < 0)
        throw std::out_of_range{"Cannot measure the given qubit"};
}

void QCircuit::reset_circuit()
{
    cached_index_ = 0;
    gate_list_.clear();
    circuit_ = af::identity(state_count(), state_count(), c32);
}

void QCircuit::generate_circuit() const
{
    if (cached_index_ != gate_list_.size())
    {
        for (std::size_t i = cached_index_; i < gate_list_.size(); ++i)
        {
            const auto& gate = *(gate_list_[i]);
            gate(const_cast<QCircuit&>(*this));
        }
        cached_index_ = gate_list_.size();
    }
}

QSimulator::QSimulator(uint32_t qubit_count, const QState& initial_state, const QNoise& noise_generator)
    :states_(qubit_count, initial_state), global_state_{fast_pow2(qubit_count), c32}, noise_{noise_generator}, qubits_{qubit_count}
{
    generate_global_state();
}

void QSimulator::generate_global_state()
{
    std::vector<af::cfloat> temp(qubit_count() * 2);
    for (uint32_t i = 0; i < qubit_count(); ++i)
    {
        temp[i * 2    ] = states_[i][0];
        temp[i * 2 + 1] = states_[i][1];
    }

    af::array states = af::array(2, qubit_count(), temp.data());

    global_state_ = states.col(0);

    //Generate the global state from the tensor product of elements
    for (uint32_t i = 1; i < qubit_count(); ++i)
        global_state_ = tensor_product(global_state_, states.col(i));
}

void QSimulator::simulate(const QCircuit& circuit)
{
    if (circuit.qubit_count() != qubits_)
        throw std::invalid_argument{"Number of qubit states and circuit input qubit states do not match"};
    circuit.generate_circuit();
    global_state_ = af::matmul(circuit.circuit(), global_state_);
}

bool QSimulator::peek_measure(uint32_t qubit) const
{
    if (qubit >= qubit_count())
        throw std::out_of_range{"Cannot measure the state of the given qubit"};

    float val = lin_dist(rd_gen);
    float prob1 = qubit_probability_true(qubit);

    return val < prob1;
}

bool QSimulator::measure(uint32_t qubit)
{
    if (qubit >= qubit_count())
        throw std::out_of_range{"Cannot measure the state of the given qubit"};

    float val = lin_dist(rd_gen);
    
    af::array probabilities = af::real(global_state_ * af::conjg(global_state_));
    af::array indices = af::iota(state_count(), 1, u32) & fast_pow2(qubit_count() - qubit - 1);

    af::array states1 = af::select(indices.as(b8), global_state_, 0.f);
    af::array states0 = af::select(indices.as(b8), 0.f, global_state_);

    float prob1 = af::sum<float>(states1 * af::conjg(states1));
    bool measurement = val < prob1;

    if (measurement)
        global_state_ = states1 / std::sqrt(prob1);
    else
        global_state_ = states0 / std::sqrt(1.f - prob1);

    return measurement;
}

uint32_t QSimulator::peek_measure_all() const
{
    // Generate randomly a value in range [0, 1]
    float val = lin_dist(rd_gen);

    af::array conj_gs = af::conjg(global_state_);

    af::array probabilities = af::real(global_state_ * conj_gs);

    //Find all the states that posses a probability greater than or equal to the random value generated
    af::array prob = af::where(af::ceil(af::accum(probabilities) - val));

    //The entry of that first value is the state measured
    uint32_t state = prob.isempty() ? 0 : prob(0).scalar<uint32_t>();

    return state;
}

uint32_t QSimulator::measure_all()
{
    uint32_t measurement = peek_measure_all();

    //Set the global state to the result of the measurement
    global_state_(af::span) = 0.f;
    global_state_(measurement) = 1.f;

    return measurement;
}


float QSimulator::qubit_probability_true(uint32_t qubit) const
{
    if (qubit >= qubit_count() || qubit < 0)
        throw std::out_of_range{"Cannot obtain probability of the given qubit"};

    float val = lin_dist(rd_gen);
    
    af::array probabilities = af::real(global_state_ * af::conjg(global_state_));
    af::array indices = af::iota(state_count(), 1, u32) & fast_pow2(qubit_count() - qubit - 1);

    af::array states1 = af::select(indices.as(b8), global_state_, 0.);

    float prob1 = af::sum<float>(states1 * af::conjg(states1));

    return prob1;
}

float QSimulator::state_probability(uint32_t state) const
{
    if (state >= state_count() || state < 0)
        throw std::out_of_range{"Cannot obtain probability of the given state"};

    af::cfloat val = global_state_(state).scalar<af::cfloat>();
    return val.real * val.real + val.imag * val.imag;
}

std::vector<uint32_t> QSimulator::profile_measure_all(uint32_t rep_count) const
{
    const auto states = state_count();
    std::vector<uint32_t> count(states);

    //Generate rep_count amount of random numbers in [0, 1] and repeat them along dim 1
    af::array rnd = af::tile(af::randu(rep_count, f32, intern_rnd_engine), 1, states);

    //Generate the cumulative probability of each state along dim 1 and repeat them along dim 0
    af::array probabilities = af::tile(af::transpose(af::accum(af::real(global_state_ * af::conjg(global_state_)))), rep_count, 1);

    //Find all the states that posses a probability greater than or equal to the random value generated as indices + 1
    auto temp = af::sum(af::ceil(probabilities - rnd), 1).as(u32);
    
    //Sort the indices of the states measured
    af::array prob = af::sort(af::join(0, temp, af::iota(states, 1, u32) + 1), 0, false);

    //Count the indices of the states measured
    af::array keys, vals;
    af::countByKey(keys, vals, prob, af::constant(1, rep_count + states));
    vals -= 1;

    //Copy it to the host
    vals.host(count.data());

    return count;
}

std::array<uint32_t , 2> QSimulator::profile_measure(uint32_t qubit, uint32_t rep_count) const
{
    if (qubit >= qubit_count())
        throw std::out_of_range{"Cannot profile measurement of the given qubit"};

    auto states = state_count();
    int32_t qubit_state = 1 << (qubit_count() - qubit - 1);

    //Generate rep_count of random floats between 0 and 1
    af::array rnd = af::randu(rep_count, f32, intern_rnd_engine);

    //Find the probability of all states
    af::array probabilities = af::real(global_state_ * af::conjg(global_state_));

    //Mark the states where the qubit is in state 1
    af::array one_indices = (af::iota(states, 1, u32) & qubit_state).as(b8);

    //Obtain the probability of the states with qubit in state 1
    af::array one_vals = af::select(one_indices, probabilities, 0LL);

    //Obtain the probability of being in state 1 and tile it with rep_count
    af::array one_sum = af::tile(af::sum(one_vals), rep_count);

    //Find all the random numbers which result in state |1> and count them
    auto temp = af::sum<uint32_t>(af::ceil(one_sum - rnd), 0);

    return { rep_count - temp , temp };
}

QCircuit& X::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto col_indices = af::iota(states, 1, s32) ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto matrix_x = af::sparse(states, states, af::constant(af::cfloat{1.f , 0.f}, states), row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_x, circuit);

    return qc;
}

/*
template<>
QCircuit& operator<<<X>(QCircuit& qc, const std::vector<X>& gates)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (gates.size() > qubits)
        throw std::invalid_argument{"Cannot add more than circuit qubit count (concurrent) gates"};

    std::array<bool, max_qubit_count> qubit_gates{};
    for (const auto& gate : gates)
    {
        auto index = gate.target_qubit;
        if (index >= qubits)
            throw std::invalid_argument{"Cannot add gate at the given qubit position"};
        if (qubit_gates[index])
            throw std::invalid_argument{"Cannot add (concurrent) gate to the same qubit multiple times"};
        qubit_gates[index] = true;
    }

    af::array identity_matrix = af::identity(2, 2, c32);
    af::array gates_matrix = qubit_gates[0] ? x_matrix : identity_matrix;

    for (uint32_t i = 1; i < qubits; ++i)
    {
        const auto& has_gate = qubit_gates[i];
        if (has_gate)
            gates_matrix = tensor_product(gates_matrix, x_matrix);
        else
            gates_matrix = tensor_product(gates_matrix, identity_matrix);
    }

    auto& circuit = qc.circuit();
    circuit = af::matmul(gates_matrix, circuit);

    return qc;
}
*/

std::string X::to_string() const
{
    return "X," + std::to_string(target_qubit) + ";";
}

QCircuit& Y::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota = af::iota(states, 1, s32);
    auto col_indices = iota ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto values = af::constant(af::cfloat{ 0.f , -1.f }, states);
    af::replace(values, (iota & (1 << target_qubit)).as(b8), af::constant(af::cfloat{ 1.f , 0.f }, states));

    auto matrix_y = af::sparse(states, states, values, row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_y, circuit);

    return qc;
}

/*
template<>
QCircuit& operator<<<Y>(QCircuit& qc, const std::vector<Y>& gates)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (gates.size() > qubits)
        throw std::invalid_argument{"Cannot add more than circuit qubit count (concurrent) gates"};

    std::array<bool, max_qubit_count> qubit_gates{};
    for (const auto& gate : gates)
    {
        auto index = gate.target_qubit;
        if (index >= qubits)
            throw std::invalid_argument{"Cannot add gate at the given qubit position"};
        if (qubit_gates[index])
            throw std::invalid_argument{"Cannot add (concurrent) gate to the same qubit multiple times"};
        qubit_gates[index] = true;
    }

    af::array identity_matrix = af::identity(2, 2, c32);
    af::array gates_matrix = qubit_gates[0] ? y_matrix : identity_matrix;

    for (uint32_t i = 1; i < qubits; ++i)
    {
        const auto& has_gate = qubit_gates[i];
        if (has_gate)
            gates_matrix = tensor_product(gates_matrix, y_matrix);
        else
            gates_matrix = tensor_product(gates_matrix, identity_matrix);
    }

    auto& circuit = qc.circuit();
    circuit = af::matmul(gates_matrix, circuit);

    return qc;
}
*/

std::string Y::to_string() const
{
    return "Y," + std::to_string(target_qubit) + ";";
}

QCircuit& Z::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota = af::iota(states, 1, s32);
    auto col_indices = iota;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto vals = af::constant(af::cfloat{ -1.f , 0.f }, states);
    af::replace(vals, (iota & mask).as(b8), af::constant(af::cfloat{ 1.f , 0.f }, states));

    auto matrix_y = af::sparse(states, states, vals, row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_y, circuit);

    return qc;
}

/*
template<>
QCircuit& operator<<<Z>(QCircuit& qc, const std::vector<Z>& gates)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (gates.size() > qubits)
        throw std::invalid_argument{"Cannot add more than circuit qubit count (concurrent) gates"};

    std::array<bool, max_qubit_count> qubit_gates{};
    for (const auto& gate : gates)
    {
        auto index = gate.target_qubit;
        if (index >= qubits)
            throw std::invalid_argument{"Cannot add gate at the given qubit position"};
        if (qubit_gates[index])
            throw std::invalid_argument{"Cannot add (concurrent) gate to the same qubit multiple times"};
        qubit_gates[index] = true;
    }

    af::array identity_matrix = af::identity(2, 2, c32);
    af::array gates_matrix = qubit_gates[0] ? z_matrix : identity_matrix;

    for (uint32_t i = 1; i < qubits; ++i)
    {
        const auto& has_gate = qubit_gates[i];
        if (has_gate)
            gates_matrix = tensor_product(gates_matrix, z_matrix);
        else
            gates_matrix = tensor_product(gates_matrix, identity_matrix);
    }

    auto& circuit = qc.circuit();
    circuit = af::matmul(gates_matrix, circuit);

    return qc;
}
*/

std::string Z::to_string() const
{
    return "Z," + std::to_string(target_qubit) + ";";
}

QCircuit& Hadamard::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

    //Find the n-qubit X matrix for the target qubit as I_L @ H @ I_R where @ is the tensor product
    af::array temp = tensor_product(left_identity, hadamard_matrix);
    temp = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit = af::matmul(temp, circuit);
    return qc;

/*
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    int mask = 1 << (qubits - target_qubit - 1);
    int masks_val[] = {
        0, mask
    };
    auto row_indices = af::iota(states + 1, 1, s32) * 2;
    auto temp = af::flat(af::tile(af::iota(states, 1, s32).T(), 2)) & af::constant(~mask, states * 2, s32);
    temp.eval();
    auto temp2 = af::tile(af::array(2, 1, masks_val), states);
    temp2.eval();

    const float sqrt2 = 0.70710678118f;
    auto values = af::constant(af::cfloat{ sqrt2 , 0.f }, states * 2);
    int minus_mask = 1 << (qubits - target_qubit) | 1;
    values.eval();
    af::replace(values, ((af::iota(states * 2, 1, s32) & minus_mask) ^ minus_mask).as(b8),
                af::constant(af::cfloat{ -sqrt2 , 0}, states * 2));

    row_indices.eval();
    column_indices.eval();
    values.eval();
    auto matrix_hadamard = af::sparse(states, states, values, row_indices, column_indices);
    matrix_hadamard.eval();
    //af_print(af::dense(matrix_hadamard));

    auto& circuit = qc.circuit();
    circuit.eval();
    circuit = af::matmul(matrix_hadamard, circuit);

    return qc;
*/
}

/*
template<>
QCircuit& operator<<<Hadamard>(QCircuit& qc, const std::vector<Hadamard>& gates)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (gates.size() > qubits)
        throw std::invalid_argument{"Cannot add more than circuit qubit count (concurrent) gates"};

    std::array<bool, max_qubit_count> qubit_gates{};
    for (const auto& gate : gates)
    {
        auto index = gate.target_qubit;
        if (index >= qubits)
            throw std::invalid_argument{"Cannot add gate at the given qubit position"};
        if (qubit_gates[index])
            throw std::invalid_argument{"Cannot add (concurrent) gate to the same qubit multiple times"};
        qubit_gates[index] = true;
    }

    af::array identity_matrix = af::identity(2, 2, c32);
    af::array gates_matrix = qubit_gates[0] ? hadamard_matrix : identity_matrix;

    for (uint32_t i = 1; i < qubits; ++i)
    {
        const auto& has_gate = qubit_gates[i];
        if (has_gate)
            gates_matrix = tensor_product(gates_matrix, hadamard_matrix);
        else
            gates_matrix = tensor_product(gates_matrix, identity_matrix);
    }

    auto& circuit = qc.circuit();
    circuit = af::matmul(gates_matrix, circuit);

    return qc;
}
*/

std::string Hadamard::to_string() const
{
    return "H," + std::to_string(target_qubit) + ";";
}

QCircuit& Phase::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const int32_t mask = 1 << (qubits - target_qubit - 1);
    auto iota = af::iota(states, 1, s32);
    auto col_indices = iota;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto vals = af::constant(af::cfloat{ std::cos(angle) , std::sin(angle) }, states);
    af::replace(vals, (iota & mask).as(b8), af::constant(af::cfloat{ 1.f , 0.f }, states));

    auto matrix_phase = af::sparse(states, states, vals, row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_phase, circuit);

    return qc;
}

/*
template<>
QCircuit& operator<<<Phase>(QCircuit& qc, const std::vector<Phase>& gates)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    if (gates.size() > qubits)
        throw std::invalid_argument{"Cannot add more than circuit qubit count (concurrent) gates"};

    std::array<bool, max_qubit_count> qubit_gates{};
    std::array<af::cfloat, max_qubit_count> qubit_angles{};
    for (const auto& gate : gates)
    {
        auto index = gate.target_qubit;
        if (index >= qubits)
            throw std::invalid_argument{"Cannot add gate at the given qubit position"};
        if (qubit_gates[index])
            throw std::invalid_argument{"Cannot add (concurrent) gate to the same qubit multiple times"};
        qubit_gates[index] = true;
        qubit_angles[index] = { std::cos(gate.angle) , std::sin(gate.angle) };
    }

    af::array identity_matrix = af::identity(2, 2, c32);
    af::array phase_matrices = af::tile(af::flat(identity_matrix), 1, max_qubit_count);
    phase_matrices.row(3) = af::array(1, max_qubit_count, qubit_angles.data());

    af::array gates_matrix = qubit_gates[0] ? af::moddims(phase_matrices.col(0), 2, 2) : identity_matrix;

    for (uint32_t i = 1; i < qubits; ++i)
    {
        const auto& has_gate = qubit_gates[i];
        if (has_gate)
            gates_matrix = tensor_product(gates_matrix, af::moddims(phase_matrices.col(i), 2, 2));
        else
            gates_matrix = tensor_product(gates_matrix, identity_matrix);
    }

    auto& circuit = qc.circuit();
    circuit = af::matmul(gates_matrix, circuit);

    return qc;
}
*/

std::string Phase::to_string() const
{
    if (angle == pi / 2)
        return "S," + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 2)
        return "S†," + std::to_string(target_qubit) + ";";
    else if (angle == pi / 4)
        return "T," + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 4)
        return "T†," + std::to_string(target_qubit) + ";";
    else
        return "Phase," + std::to_string(target_qubit) + ";";
}


QCircuit& Swap::operator()(QCircuit& qc) const
{
    const auto states = qc.state_count();
    const auto qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{"Cannot use the swap gate on the same target qubits"};

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int32_t maskA = 1 << (qubits - 1 - target_qubit_A);
    int32_t maskB = 1 << (qubits - 1 - target_qubit_B);

    // Generate column indices (swap the bits in the states of qA and qB)
    auto col_indices = af::iota(states, 1, s32);
    auto temp = (col_indices >> posA) ^ (col_indices >> posB);
    temp = temp & 1;
    col_indices = col_indices ^ ((temp << posA) | (temp << posB));

    auto row_indices = af::iota(states + 1, 1, s32);

    auto swap_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.f }, states), row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(swap_matrix, circuit);

    return qc;
}

std::string Swap::to_string() const
{
    return "Swap," + std::to_string(target_qubit_A) + "," + std::to_string(target_qubit_B) + ";";
}

QCircuit& Control_X::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8), iota ^ target_mask);

    auto cx_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.f , 0.f }, states), row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(cx_matrix, circuit);

    return qc;
}

std::string Control_X::to_string() const
{
    return "CX,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& Control_Y::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);
    int32_t mask = control_mask | target_mask;

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8), iota ^ target_mask);

    auto values = af::constant(af::cfloat{ 1.f , 0.f }, states);
    auto op_indices = iota & mask;
    af::replace(values, op_indices != mask, af::constant(af::cfloat{ 0.f , 1.f }, states));
    af::replace(values, op_indices != control_mask, af::constant(af::cfloat{ 0.f , -1.f }, states));

    auto cx_matrix = af::sparse(states, states, values, row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(cx_matrix, circuit);

    return qc;
}

std::string Control_Y::to_string() const
{
    return "CY,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& Control_Z::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);
    int32_t mask = control_mask | target_mask;

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    auto values = af::constant(af::cfloat{ 1.f , 0.f }, states);
    af::replace(values, (iota & mask) != mask, af::constant(af::cfloat{ -1.f , 0.f }, states));
    //af::replace(values, ((iota & mask) ^ mask).as(b8), af::constant(af::cfloat{ -1.f , 0.f}, states));

    auto cz_matrix = af::sparse(states, states, values, row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(cz_matrix, circuit);

    return qc;
}

std::string Control_Z::to_string() const
{
    return "CZ,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& Control_Phase::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);
    int32_t mask = control_mask | target_mask;

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;

    auto values = af::constant(af::cfloat{ 1.f , 0.f }, states);
    af::replace(values, (iota & mask) != mask, af::constant(af::cfloat{ std::cos(angle) , std::sin(angle) }, states));

    auto cphase_matrix = af::sparse(states, states, values, row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(cphase_matrix, circuit);

    return qc;
}

std::string Control_Phase::to_string() const
{
    if (angle == pi / 2)
        return "CS,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 2)
        return "CS†,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
    else if (angle == pi / 4)
        return "CT,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
    else if (angle == -pi / 4)
        return "CT†,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
    else
        return "CPhase,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

static QCircuit& cswap_gate_af_cpu(QCircuit& qc, uint32_t control_qubit, uint32_t target_qubit_A, uint32_t target_qubit_B)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    std::vector<int> cols(states), rows(states + 1);

    int32_t control_mask = 1 << (qubits - 1 - control_qubit);
    int32_t target_maskA = 1 << (qubits - 1 - target_qubit_A);
    int32_t target_maskB = 1 << (qubits - 1 - target_qubit_B);
    int32_t target_mask = target_maskA | target_maskB;

    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        rows[i + 1] = i + 1;
        if ((i & control_mask) && ((i & target_mask) == target_maskA || (i & target_mask) == target_maskB))
        {
            cols[i] = (i ^ target_maskA) ^ target_maskB;
        }
        else
            cols[i] = i;
    }

    af::array cswap_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                            af::array(states + 1, rows.data()), af::array(states, cols.data()));

    auto& circuit = qc.circuit();
    circuit = af::matmul(cswap_matrix, circuit);

    return qc;
}

static QCircuit& cswap_gate_af_opencl(QCircuit& qc, uint32_t control_qubit, uint32_t target_qubit_A, uint32_t target_qubit_B)
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int32_t control = 1 << (qubits - control_qubit - 1);
    int32_t targetA = 1 << posA;
    int32_t targetB = 1 << posB;

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);

    // Generate column indices that should be swapped (swap the bits in the states of qA and qB)
    auto temp = (iota >> posA) ^ (iota >> posB);
    temp = temp & 1;
    auto replace_vals = iota ^ ((temp << posA) | (temp << posB));

    auto column_indices = iota;
    af::replace(column_indices, (iota & control) != control, replace_vals);

    auto matrix_cswap = af::sparse(states, states, af::constant(af::cfloat{ 1.f , 0.f }, states), row_indices, column_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_cswap, circuit);

    return qc;
}

QCircuit& Control_Swap::operator()(QCircuit& qc) const
{
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
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};
    if (target_qubit_A == target_qubit_B)
        throw std::invalid_argument{"Cannot use the swap gate on the same target qubits"};

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int32_t control = 1 << (qubits - control_qubit - 1);
    int32_t targetA = 1 << posA;
    int32_t targetB = 1 << posB;

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);

    // Generate column indices that should be swapped (swap the bits in the states of qA and qB)
    auto temp = (iota >> posA) ^ (iota >> posB);
    temp = temp & 1;
    auto replace_vals = iota ^ ((temp << posA) | (temp << posB));

    auto column_indices = iota;
    af::replace(column_indices, (iota & control) != control, replace_vals);

    auto matrix_cswap = af::sparse(states, states, af::constant(af::cfloat{ 1.f , 0.f }, states), row_indices, column_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_cswap, circuit);

    return qc;
}

std::string Control_Swap::to_string() const
{
    return "CSwap" + std::to_string(control_qubit) + "," + std::to_string(target_qubit_A) + 
            std::to_string(target_qubit_B) + ";";
}

QCircuit& Control_Hadamard::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit == target_qubit)
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    const int32_t value_mask = (1 << (qubits - target_qubit)) | 1;
    const int32_t control_index_mask = 1 << (qubits - control_qubit);
    const int32_t index_mask = ~(1 << (qubits - target_qubit - 1));

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32) * 2;
    //auto column_indices = af::resize(iota, states * 2, 1, AF_INTERP_LOWER) & index_mask;
    auto resized_iota = af::flat(af::tile(iota.T(), 2));
    auto column_indices = resized_iota & index_mask;

    const int32_t off[] = {
        0, 1 << (qubits - target_qubit - 1)
    };

    auto offset = af::tile(af::array(2, off), states);
    column_indices += offset;

    auto iota2 = af::iota(states * 2, 1, s32);
    auto control_indices = (iota2 & control_index_mask).as(b8);
    //af::replace(column_indices, control_indices, af::resize(iota, states * 2, 1, AF_INTERP_NEAREST));
    auto new_col_indices = af::shift(resized_iota, -1);
    new_col_indices(af::end) = states;
    af::replace(column_indices, control_indices, new_col_indices);

    const float sqrt2 = 0.70710678118f;
    const int32_t val_mask = (1 << (qubits - target_qubit)) | 1;

    const af::cfloat identity_vals[] = {
        { 1.f , 0.f } , { 0.f , 0.f }
    };

    auto values = af::tile(af::array(2, identity_vals), states);
    auto value_indices = ((iota2 & value_mask) ^ value_mask).as(b8);

    af::replace(values, !(control_indices &&  value_indices), af::constant(af::cfloat{ sqrt2, 0.f }, states * 2));
    af::replace(values, !(control_indices && !value_indices), af::constant(af::cfloat{ -sqrt2, 0.f }, states * 2));

    auto matrix_chadamard = af::sparse(states, states, values, row_indices, column_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_chadamard, circuit);

    return qc;
}

std::string Control_Hadamard::to_string() const
{
    return "CH," + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& CControl_Not::operator()(QCircuit& qc) const
{
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
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;
    af::replace(column_indices, ((iota ^ control_mask) & control_mask).as(b8), iota ^ target_mask);

    auto ccnot_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.f , 0.f }, states), row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(ccnot_matrix, circuit);

    return qc;
}

std::string CControl_Not::to_string() const
{
    return "CCX,2,1:" + std::to_string(control_qubit_A) + "," + std::to_string(control_qubit_B) + 
            "," + std::to_string(target_qubit) + ";";
}

QCircuit& Or::operator()(QCircuit& qc) const
{
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
        throw std::invalid_argument{"Control qubit cannot be the same as the target qubit"};

    int32_t control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
    int32_t target_mask = 1 << (qubits - 1 - target_qubit);

    auto iota = af::iota(states, 1, s32);
    auto row_indices = af::iota(states + 1, 1, s32);
    auto column_indices = iota;
    af::replace(column_indices, ((iota & control_mask) == 0).as(b8), iota ^ target_mask);

    auto ccnot_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.f , 0.f }, states), row_indices, column_indices);

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(ccnot_matrix, circuit);

    return qc;
}

std::string Or::to_string() const
{
    auto top = control_qubit_A < control_qubit_B ? control_qubit_A : control_qubit_B;
    auto bottom = control_qubit_A < control_qubit_B ? control_qubit_B : control_qubit_A;
    if (target_qubit < top)
        top = target_qubit;
    if (target_qubit > bottom)
        bottom = target_qubit;

    std::string qubits;
    for (uint32_t i = top; i < bottom; i++)
        qubits.append(std::to_string(i) + ",");
    qubits.append(std::to_string(bottom) + ";");

    return "Or,0," + std::to_string(bottom - top + 1) + ":" + qubits;
}

static std::string update_circuit_representation(const std::string& circuit_string, uint32_t offset);
static std::string update_ctrl_circuit_representation(const std::string& circuit_string, uint32_t control, uint32_t target);

CircuitGate::CircuitGate(const QCircuit& circuit_, uint32_t target_qubit_begin_, std::string name)
    : internal_circuit(circuit_), representation{},
        qubit_count{circuit_.qubit_count()}, target_qubit_begin{target_qubit_begin_}
{
    if (name == "")
    {
        representation = update_circuit_representation(circuit_.representation(), target_qubit_begin_);
    }
    else
    {
        //Naming of the circuit must not contain parsing tokens
        if (name.find_first_of(",;:") != std::string::npos)
            throw std::invalid_argument{"Name cannot contain commas, colons, nor semicolons"};

        //Set circuit to have no ctrl qubits
        representation.append(name).append(",0,").append(std::to_string(circuit_.qubit_count())).append(":");

        //Offset the target qubits by the target_qubit_begin
        for (int i = 0; i < circuit_.qubit_count() - 1; ++i)
            representation.append(std::to_string(i + target_qubit_begin_)).append(",");
        representation.append(std::to_string(circuit_.qubit_count() - 1 + target_qubit_begin_));
        representation.append(";");
    }
}


QCircuit& CircuitGate::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    const auto state_count = fast_pow2(qubit_count);
    if (target_qubit_begin >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    const uint32_t circuit_qubits = qc.qubit_count();
    const uint32_t circuit_states = qc.state_count();
    const uint32_t gate_qubits = qubit_count;
    const uint32_t gate_states = 1 << gate_qubits;
    const uint32_t gate_qubit_begin = target_qubit_begin;
    af::array gate_matrix = af::identity(circuit_states, circuit_states, c32);

    internal_circuit.generate_circuit();
    const af::array& gate = internal_circuit.circuit();

    uint32_t rem_count = 1 << (circuit_qubits - gate_qubits);
    
    auto len = gate_states * gate_states * rem_count;
    auto m = af::tile(af::flat(af::tile(af::iota(gate_states, 1, s32).T(), gate_states)), rem_count);
    auto n = af::iota(gate_states, gate_states * rem_count, s32);
    auto ind = af::flat(af::tile(af::iota(rem_count, 1, s32).T(), gate_states * gate_states));
    auto ii = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, n, ind, len);
    auto jj = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, m, ind, len);

    af::array gate_values = gate(n * gate_states + m);

    jj.eval();
    //af_print(ii.as(f32) * af::constant(16, len, f32)+ jj.as(f32));
    gate_matrix(ii * circuit_states + jj) = gate_values;

    /*
    auto len = gate_states * gate_states;
    auto m = af::iota(gate_states, gate_states, s32);
    auto n = af::resize(af::iota(gate_states, 1, s32), len, 1, AF_INTERP_LOWER);
    af::array gate_values = gate(m * gate_states + n);
    for (int i = 0; i < rem_count; ++i)
    {
        auto ind = af::constant(i, len, s32);
        auto ii = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, m, ind, len);
        auto jj = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, n, ind, len);
        gate_matrix(ii * circuit_states + jj) = gate_values;
    }*/

    auto& circuit = qc.circuit();
    circuit = af::matmul(gate_matrix, circuit);

    return qc;
}

ControlCircuitGate::ControlCircuitGate(const QCircuit& circuit_, uint32_t control_qubit_, uint32_t target_qubit_begin_, std::string name)
        : internal_circuit(circuit_), representation{},
          qubit_count{circuit_.qubit_count()}, control_qubit{control_qubit_}, target_qubit_begin{target_qubit_begin_}
{
    if (name == "")
    {
        representation = update_ctrl_circuit_representation(circuit_.representation(), control_qubit, target_qubit_begin_);
    }
    else
    {
        //Naming of the circuit must not contain parsing tokens
        if (name.find_first_of(",;:") != std::string::npos)
            throw std::invalid_argument{"Name cannot contain commas, colons, nor semicolons"};

        //Set the circuit to contain 1 ctrl qubit and append it to the front of the qubit list
        representation.append(name).append(",1,").append(std::to_string(circuit_.qubit_count())).append(":");
        representation.append(std::to_string(control_qubit_)).append(",");

        //Offset all the target qubits by target_qubit_begin
        for (int i = 0; i < circuit_.qubit_count() - 1; ++i)
            representation.append(std::to_string(i + target_qubit_begin_)).append(",");
        representation.append(std::to_string(circuit_.qubit_count() - 1 + target_qubit_begin_));
        representation.append(";");
    }
}

QCircuit& ControlCircuitGate::operator()(QCircuit& qc) const
{
    const auto qubits = qc.qubit_count();
    const auto states = qc.state_count();
    const auto state_count = fast_pow2(qubit_count);

    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Gate must fit inside the circuit qubit count"};
    if (target_qubit_begin <= control_qubit && control_qubit < target_qubit_begin + qubit_count)
        throw std::out_of_range{"Control qubit cannot be one of the target qubits of the gate"};

    const uint32_t circuit_qubits = qc.qubit_count();
    const uint32_t circuit_states = qc.state_count();
    const uint32_t gate_qubits = qubit_count;
    const uint32_t gate_states = 1 << gate_qubits;
    const uint32_t gate_qubit_begin = target_qubit_begin;
    af::array gate_matrix = af::identity(circuit_states, circuit_states, c32);

    internal_circuit.generate_circuit();
    const af::array& gate = internal_circuit.circuit();

    uint32_t rem_count = 1 << (circuit_qubits - gate_qubits - 1);
    
    auto len = gate_states * gate_states * rem_count;
    auto m = af::tile(af::flat(af::tile(af::iota(gate_states, 1, s32).T(), gate_states)), rem_count);
    auto n = af::iota(gate_states, gate_states * rem_count, s32);
    auto ind = af::flat(af::tile(af::iota(rem_count, 1, s32).T(), gate_states * gate_states));
    auto ii = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, n, ind, len);
    auto jj = gen_index(gate_qubit_begin, gate_qubits, circuit_qubits, m, ind, len);

    auto ones = af::constant(1, len, s32);
    ii = insert_bits(ii, ones, af::constant(qubits - control_qubit - 1, len, s32), ones, len);
    jj = insert_bits(jj, ones, af::constant(qubits - control_qubit - 1, len, s32), ones, len);

    af::array gate_values = gate(n * gate_states + m);

    jj.eval();
    gate_matrix(ii * circuit_states + jj) = gate_values;

    auto& circuit = qc.circuit();
    circuit = af::matmul(gate_matrix, circuit);

    return qc;
}

void initialize(int argc, char** argv)
{
    int device = (argc > 1) ? ::std::stoi(argv[1]) : 0;
    af::setDevice(device);
    af::info();

    x_matrix = af::array(2, 2, x_mat);

    y_matrix = af::array(2, 2, y_mat);

    z_matrix = af::array(2, 2, z_mat);

    hadamard_matrix = af::array(2, 2, h_mat);

    std::random_device rnd_device;
    intern_rnd_engine = af::randomEngine(AF_RANDOM_ENGINE_THREEFRY, rnd_device());
}

QState X_op(const QState& state)
{
    af::cfloat temp[2];
    temp[0] = (x_mat[0] * state[0]) + (x_mat[1] * state[1]);
    temp[1] = (x_mat[2] * state[0]) + (x_mat[3] * state[1]);
    
    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}
    
QState Y_op(const QState& state)
{
    af::cfloat temp[2];
    temp[0] = (y_mat[0] * state[0]) + (y_mat[1] * state[1]);
    temp[1] = (y_mat[2] * state[0]) + (y_mat[3] * state[1]);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState Z_op(const QState& state)
{
    af::cfloat temp[2];
    temp[0] = (z_mat[0] * state[0]) + (z_mat[1] * state[1]);
    temp[1] = (z_mat[2] * state[0]) + (z_mat[3] * state[1]);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState RotateX_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * std::cos(angle / 2.f) + state[1] * af::cfloat{ 0.f , -std::sin(angle / 2.f) };
    temp[1] = state[0] * af::cfloat{ 0.f , -std::sin(angle / 2.f) } + state[1] * std::cos(angle / 2.f);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState RotateY_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * std::cos(angle / 2.f) + state[1] * -std::sin(angle / 2.f);
    temp[1] = state[0] * std::sin(angle / 2.f) + state[1] *  std::cos(angle / 2.f);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState RotateZ_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * af::cfloat{ std::cos(angle / 2.f) , -std::sin(angle / 2.f)} + state[1] * 0.f;
    temp[1] = state[0] * 0.f + state[1] * af::cfloat{ std::cos(angle / 2.f) , std::sin(angle / 2.f)};

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState Hadamard_op(const QState& state)
{
    af::cfloat temp[2];
    temp[0] = h_mat[0] * state[0] + h_mat[1] * state[1];
    temp[1] = h_mat[2] * state[0] + h_mat[3] * state[1];

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState Phase_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * 1.f  + state[1] * 0.f;
    temp[1] = state[0] * 0.f  + state[1] * af::cfloat{std::cos(angle) , std::sin(angle)};

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

std::string update_circuit_representation(const std::string& circuit_string, uint32_t offset)
{
    std::string out;
    std::size_t current_begin = 0;
    auto current_end = circuit_string.find(";", current_begin);

    while (current_end != std::string::npos)
    {
        auto name_end = circuit_string.find(",", current_begin);
        auto name_str = circuit_string.substr(current_begin, name_end - current_begin);
        if (name_end > current_end)
        {
            out.append(circuit_string.substr(current_begin, current_end - current_begin));
            out.append(";");
            current_begin = current_end + 1;
            current_end = circuit_string.find(";", current_begin);
            continue;
        }

        out.append(name_str).append(",");

        auto colon_pos = circuit_string.find(":", current_begin);
        auto qubit_pos_begin = name_end + 1;
        if (colon_pos < current_end)
        {
            out.append(circuit_string.substr(name_end + 1, colon_pos - name_end));
            qubit_pos_begin = colon_pos + 1;
        }

        auto qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        qubit_pos_end = qubit_pos_end > current_end ? current_end : qubit_pos_end;

        //Add the offset to all target qubits
        while (qubit_pos_end < current_end)
        {
            out.append(std::to_string(std::stoi(circuit_string.substr(qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) + offset));
            out.append(",");
            qubit_pos_begin = qubit_pos_end + 1;
            qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        }
        qubit_pos_end = qubit_pos_end > current_end ? current_end : qubit_pos_end;
        out.append(std::to_string(std::stoi(circuit_string.substr(qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) + offset));

        current_begin = current_end + 1;
        current_end = circuit_string.find(";", current_begin);
        out.append(";");
    }

    return out;
}

std::string update_ctrl_circuit_representation(const std::string& circuit_string, uint32_t control, uint32_t target)
{
    std::string out;
    std::size_t current_begin = 0;
    auto current_end = circuit_string.find(";", current_begin);

    while (current_end != std::string::npos)
    {
        auto name_end = circuit_string.find(",", current_begin);
        auto name_str = circuit_string.substr(current_begin, name_end - current_begin);
        if (name_end > current_end)
        {
            out.append(circuit_string.substr(current_begin, current_end - current_begin));
            out.append(";");
            current_begin = current_end + 1;
            current_end = circuit_string.find(";", current_begin);
            continue;
        }
        out.append("C" + name_str).append(",");

        //If matrix is fundamental add token for custom gate and set ctrl qubit count to 1
        //Else add 1 to the ctrl qubit count to the gate
        auto colon_pos = circuit_string.find(":", current_begin);
        auto qubit_pos_begin = name_end + 1;
        if (colon_pos < current_end)
        {
            auto temp = circuit_string.find(",", name_end + 1);
            out.append(std::to_string(std::stoi(circuit_string.substr(name_end + 1, temp - name_end - 1)) + 1));
            out.append(circuit_string.substr(temp, colon_pos - temp + 1));
            qubit_pos_begin = colon_pos + 1;
        }
        else
        {
            if (name_str == "X" || name_str == "Y" || name_str == "Z" || name_str == "Phase")
                out.append("1,1:");
            else if (name_str == "Swap")
                out.append("1,2:");
            // else if (name_str == "CX" || name_str == "CY" || name_str == "CZ" || name_str == "CPhase")
            //     out.append("2,1:");
        }

        out.append(std::to_string(control)).append(",");

        auto qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        qubit_pos_end = qubit_pos_end > current_end ? current_end : qubit_pos_end;

        //Add the offset to all target qubits
        while (qubit_pos_end < current_end)
        {
            out.append(std::to_string(std::stoi(circuit_string.substr(qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) + target));
            out.append(",");
            qubit_pos_begin = qubit_pos_end + 1;
            qubit_pos_end = circuit_string.find(",", qubit_pos_begin);
        }
        qubit_pos_end = qubit_pos_end > current_end ? current_end : qubit_pos_end;
        out.append(std::to_string(std::stoi(circuit_string.substr(qubit_pos_begin, qubit_pos_end - qubit_pos_begin)) + target));

        current_begin = current_end + 1;
        current_end = circuit_string.find(";", current_begin);
        out.append(";");
    }

    return out;
}

}