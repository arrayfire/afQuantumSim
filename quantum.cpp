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
#include <iostream>
#include <iomanip>

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
    Y = [  0 , i ]
        [ -i , 0 ] 
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

namespace aqs
{

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
    bool measurement = static_cast<int>(peek_measure());

    state_[ measurement] = 1.f;
    state_[!measurement] = 0.f;

    return measurement;
}

std::array<uint32_t, 2> QState::profile_measure(uint32_t rep_count) const
{
    uint32_t t = 0;
    for (int i = 0; i < rep_count; ++i)
        t += static_cast<int>(peek_measure());

    return { rep_count - t, t };
}

void QState::force_normalize()
{
    float mag2 = state_[0].real * state_[0].real + state_[0].imag * state_[0].imag +
                 state_[1].real * state_[1].real + state_[1].imag * state_[1].imag;
    
    if (mag2 == 0.f)
        throw std::invalid_argument{"Cannot normalize a null state"};
    
    float mag = sqrtf(mag2);
    
    state_[0] = state_[0] / mag;
    state_[1] = state_[1] / mag;
}

QCircuit::QCircuit(uint32_t qubit_count)
    :circuit_{af::identity(fast_pow2(qubit_count), fast_pow2(qubit_count), c32)}, representation_{}, qubits_{qubit_count}
{
    if (qubit_count < 1)
        throw std::invalid_argument{"Circuit must contain at least 1 qubit"};
    if (qubit_count > 30)
        throw std::invalid_argument{"Maximum qubit count supported is 30"};
}
    
void QCircuit::Global_Measure()
{
}

void QCircuit::Measure(uint32_t qubit)
{
    const int count = qubit_count();
    if (qubit >= count || qubit < 0)
        throw std::out_of_range{"Cannot measure the given qubit"};
}

void QCircuit::reset_circuit()
{
    circuit_ = af::identity(state_count(), state_count(), c32);
}

QSimulator::QSimulator(uint32_t qubit_count, const QState& initial_state, const QNoise& noise_generator)
    :states_(qubit_count, initial_state), global_state_{fast_pow2(qubit_count), c32}, noise_{noise_generator}, qubits_{qubit_count}
{
    generate_global_state();
}

void QSimulator::generate_global_state()
{
    std::vector<af::cfloat> temp(qubit_count() * 2);
    for (int i = 0; i < qubit_count(); ++i)
    {
        temp[i * 2    ] = states_[i][0];
        temp[i * 2 + 1] = states_[i][1];
    }

    af::array states = af::array(2, qubit_count(), temp.data());

    global_state_ = states.col(0);

    //Generate the global state from the tensor product of elements
    for (int i = 1; i < qubit_count(); ++i)
        global_state_ = tensor_product(global_state_, states.col(i));
}

void QSimulator::simulate(const QCircuit& circuit)
{
    if (circuit.qubit_count() != qubits_)
        throw std::invalid_argument{"Number of qubit states and circuit input qubit states do not match"};
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
        global_state_ = states1 / sqrtf(prob1);
    else
        global_state_ = states0 / sqrt(1.f - prob1);

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
    
    int qubit_one = fast_log2(measurement);

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
    const int states = state_count();
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

    int states = state_count();
    int qubit_state = 1 << (qubit_count() - qubit - 1);

    //Generate rep_count of random floats between 0 and 1
    af::array rnd = af::randu(rep_count, f32, intern_rnd_engine);

    //Find the probability of all states
    af::array probabilities = af::real(global_state_ * af::conjg(global_state_));

    //Mark the states where the qubit is in state 1
    af::array one_indices = (af::iota(states, 1, u32) & qubit_state).as(b8);

    //Obtain the probability of the states with qubit in state 1
    af::array one_vals = af::select(one_indices, probabilities, 0);

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

    const int mask = 1 << (qubits - target_qubit - 1);
    auto col_indices = af::iota(states, 1, s32) ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto matrix_x = af::sparse(states, states, af::constant(af::cfloat{1.f , 0.f}, states), row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_x, circuit);

    return qc;
}

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

    const int mask = 1 << (qubits - target_qubit - 1);
    auto iota = af::iota(states, 1, s32);
    auto col_indices = iota ^ mask;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto vals = af::constant(af::cfloat{ 0.f , -1.f }, states);
    af::replace(vals, (iota & (1 << target_qubit)).as(b8), af::constant(af::cfloat{ 1.f , 0.f }, states));

    auto matrix_y = af::sparse(states, states, vals, row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_y, circuit);

    return qc;
}

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

    const int mask = 1 << (qubits - target_qubit - 1);
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

std::string Z::to_string() const
{
    return "Z," + std::to_string(target_qubit) + ";";
}

QCircuit& Hadamard::operator()(QCircuit& qc) const
{
    const int qubits = qc.qubit_count();
    if (target_qubit >= qubits || target_qubit < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
    af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

    //Find the n-qubit Hadamard matrix for the target qubit as I_L @ Hadamard @ I_R where @ is the tensor product
    af::array temp = tensor_product(left_identity, hadamard_matrix);
    temp = tensor_product(temp, right_identity);

    auto& circuit = qc.circuit();
    circuit = af::matmul(temp, circuit);

    return qc;
}

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

    const int mask = 1 << (qubits - target_qubit - 1);
    auto iota = af::iota(states, 1, s32);
    auto col_indices = iota;
    auto row_indices = af::iota(states + 1, 1, s32);

    auto vals = af::constant(af::cfloat{ std::cos(angle) , std::sin(angle) }, states);
    af::replace(vals, (iota & mask).as(b8), af::constant(af::cfloat{ 1.f , 0.f }, states));

    auto matrix_y = af::sparse(states, states, vals, row_indices, col_indices);

    auto& circuit = qc.circuit();
    circuit = af::matmul(matrix_y, circuit);

    return qc;
}

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
    if (target_qubit_A < 0 || target_qubit_A >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B < 0 || target_qubit_B >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    int32_t posA = qubits - 1 - target_qubit_A;
    int32_t posB = qubits - 1 - target_qubit_B;

    int maskA = 1 << (qubits - 1 - target_qubit_A);
    int maskB = 1 << (qubits - 1 - target_qubit_B);

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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit < 0 || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit < 0 || target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);

    int control_mask = 1 << (qubits - 1 - control_qubit);
    int target_mask = 1 << (qubits - 1 - target_qubit);
    int mask = control_mask | target_mask;

    //Generate the sparse matrix positions
    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        rows[i + 1] = i + 1;
        cols[i] = (i & control_mask) ? (i ^ target_mask) : i;
    }

    //All entries are (1.f, 0.0f), generate the sparse matrix
    af::array cx_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f}, states), af::array(states + 1, rows.data()), af::array(states, cols.data()));

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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit < 0 || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit < 0 || target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);
    std::vector<af::cfloat> vals(states);

    int control_mask = 1 << (qubits - 1 - control_qubit);
    int target_mask = 1 << (qubits - 1 - target_qubit);
    int mask = control_mask | target_mask;

    // Generate the sparse matrix entry positions and values
    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        vals[i] = (i & control_mask) ? (i & target_mask ? af::cfloat{0.f, 1.f} : af::cfloat{0.f, -1.f}) : 1.0f;
        rows[i + 1] = i + 1;
        cols[i] = (i & control_mask) ? (i ^ target_mask) : i;
    }

    //Generate operation matrix
    af::array cy_matrix = af::sparse(states, states, af::array(states, vals.data()), af::array(states + 1, rows.data()),
                                       af::array(states, cols.data()));

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(cy_matrix, circuit);

    return qc;
}

std::string Control_Y::to_string() const
{
    return "CY,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& Control_Z::operator()(QCircuit& qc) const
{
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit < 0 || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit < 0 || target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);
    std::vector<af::cfloat> vals(states);

    int control_mask = 1 << (qubits - 1 - control_qubit);
    int target_mask = 1 << (qubits - 1 - target_qubit);
    int mask = control_mask | target_mask;

    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        vals[i] = (i & mask) != mask ? 1.0f : -1.0f;
        rows[i + 1] = i + 1;
        cols[i] = i;
    }

    af::array cz_matrix = af::sparse(states, states, af::array(states, vals.data()),
                                       af::array(states + 1, rows.data()), af::array(states, cols.data()));

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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit < 0 || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit < 0 || target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);
    std::vector<af::cfloat> vals(states);

    int control_mask = 1 << (qubits - 1 - control_qubit);
    int target_mask = 1 << (qubits - 1 - target_qubit);
    int mask = control_mask | target_mask;

    //Generate the sparse entries locations and values
    af::cfloat rot = {cosf(angle), sinf(angle)};
    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        vals[i] = ((i & mask) == mask)? rot : 1.0f;
        rows[i + 1] = i + 1;
        cols[i] = i;
    }

    //Generate the operation matrix
    af::array cphase_matrix = af::sparse(states, states, af::array(states, vals.data()),
                                       af::array(states + 1, rows.data()), af::array(states, cols.data()));

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

QCircuit& Control_Swap::operator()(QCircuit& qc) const
{
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (qubits < 3)
        throw std::domain_error{"Gate not supported for given circuit"};
    if (target_qubit_A >= qubits || target_qubit_A < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_B >= qubits || target_qubit_B < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit >= qubits || control_qubit < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);

    int control_mask = 1 << (qubits - 1 - control_qubit);
    int target_maskA = 1 << (qubits - 1 - target_qubit_A);
    int target_maskB = 1 << (qubits - 1 - target_qubit_B);
    int target_mask = target_maskA | target_maskB;

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

std::string Control_Swap::to_string() const
{
    return "CSwap" + std::to_string(control_qubit) + "," + std::to_string(target_qubit_A) + 
            std::to_string(target_qubit_B) + ";";
}

QCircuit& Control_Hadamard::operator()(QCircuit& qc) const
{
    const int qubits = qc.qubit_count();
    if (qubits < 2)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit < 0 || control_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit < 0 || target_qubit >= qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    QCircuit h(1);
    h << Hadamard(0);

    qc << ControlCircuitGate(h, control_qubit, target_qubit);

    return qc;
}

std::string Control_Hadamard::to_string() const
{
    return "CH," + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";";
}

QCircuit& CControl_Not::operator()(QCircuit& qc) const
{
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();

    if (qubits < 3 || control_qubit_A == target_qubit || control_qubit_B == target_qubit)
        throw std::domain_error{"Gate not supported for given simulation"};
    if (control_qubit_A >= qubits || control_qubit_A < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits || control_qubit_B < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits || target_qubit < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);

    int control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
    int target_mask = 1 << (qubits - 1 - target_qubit);

    //Generate the matrix entry indices
    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        rows[i + 1] = i + 1;
        cols[i] = (i & control_mask) == control_mask ? (i ^ target_mask) : i;
    }

    //All entries are 1, generate the operation matrix 
    af::array ccnot_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                        af::array(states + 1, rows.data()), af::array(states, cols.data()));

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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    if (control_qubit_A >= qubits || control_qubit_A < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (control_qubit_B >= qubits || control_qubit_B < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit >= qubits || target_qubit < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    std::vector<int> cols(states), rows(states + 1);

    int control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
    int target_mask = 1 << (qubits - 1 - target_qubit);

    //Generate the sparse matrix entries
    rows[0] = 0;
    for (int i = 0; i < states; ++i)
    {
        rows[i + 1] = i + 1;
        cols[i] = (i & control_mask) ? i ^ target_mask : i;
    }

    //Generate the operation matrix
    af::array or_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                     af::array(states + 1, rows.data()), af::array(states, cols.data()));

    //Update the circuit matrix by matrix multiplication
    auto& circuit = qc.circuit();
    circuit = af::matmul(or_matrix, circuit);

    return qc;
}

std::string Or::to_string() const
{
    int top = control_qubit_A < control_qubit_B ? control_qubit_A : control_qubit_B;
    int bottom = control_qubit_A < control_qubit_B ? control_qubit_B : control_qubit_A;
    if (target_qubit < top)
        top = target_qubit;
    if (target_qubit > bottom)
        bottom = target_qubit;

    std::string qubits;
    for (int i = top; i < bottom; i++)
        qubits.append(std::to_string(i) + ",");
    qubits.append(std::to_string(bottom) + ";");

    return "Or,0," + std::to_string(bottom - top + 1) + ":" + qubits;
}

static std::string update_circuit_representation(const std::string& circuit_string, int offset);
static std::string update_ctrl_circuit_representation(const std::string& circuit_string, int control, int target);

CircuitGate::CircuitGate(const QCircuit& circuit_, uint32_t target_qubit_begin_, std::string name)
    : internal_circuit{circuit_.circuit()}, representation{},
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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    const int state_count = fast_pow2(qubit_count);
    if (target_qubit_begin >= qubits || target_qubit_begin < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    //Store circuit array in the host
    af::dim4 dims = internal_circuit.dims();
    std::vector<af::cfloat> vals(dims[0] * dims[1]);
    internal_circuit.host(vals.data());

    // Inserts a given bit to the given value at the given position
    // e.g. insert_bit(5, 1, 1) 5 => b101 = b1011 = 11
    //                                 ^ 
    //                                 1 
    auto insert_bit = [](int value, bool bit, int position) {
        int temp = value & ~ ((1 << position) - 1);
        temp <<= 1;
        temp &= ~(1 << position);
        temp |= (bit ? 1 : 0) << position;
        temp |= value & ((1 << position) - 1);
        return temp;
    };

    //Generates the index of the position of the element in the larget operation circuit matrix from the target qubits of the circuit,
    //and the current index being probed of the input circuit matrix
    auto gen_index = [insert_bit](int targets, int target_count, int target_value, int qubit_count, int ii) {
        int index = ii;
        for (int i = 0; i < target_count; ++i)
            index = insert_bit(index, target_value & (1 << (target_count - 1 - i)), qubit_count - target_count - targets);
        return index;
    };

    //Create a temporary storage for the operation matrix
    int rem_count = fast_pow2(qubits - qubit_count);
    std::vector<af::cfloat> out_temp;
    out_temp.resize(states * states);

    //Generate a identity matrix
    for (int i = 0; i < states; ++i)
        out_temp[i * states + i] = 1.0f;

    //Fill all the entries with the correct values
    for (int i = 0; i < rem_count; ++i)
    {
        for (int m = 0; m < state_count; ++m)
        {
            int ii = gen_index(target_qubit_begin, qubit_count, m, qubits, i);
            for (int n = 0; n < state_count; ++n)
            {
                int jj = gen_index(target_qubit_begin, qubit_count, n, qubits, i);

                //Copy the entries from the input circuit to the correct positions
                out_temp[ii * states + jj] = vals[m * state_count + n];
            }
        }
    }

    //Generate the operation matrix
    af::array gate_mat = af::array(states, states, out_temp.data());

    //Update the circuit matrix
    auto& circuit = qc.circuit();
    circuit = af::matmul(gate_mat, circuit);

    return qc;
}

ControlCircuitGate::ControlCircuitGate(const QCircuit& circuit_, uint32_t control_qubit_, uint32_t target_qubit_begin_, std::string name)
        : internal_circuit{circuit_.circuit()}, representation{},
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
    const int qubits = qc.qubit_count();
    const int states = qc.state_count();
    const int state_count = fast_pow2(qubit_count);

    if (control_qubit >= qubits || control_qubit < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin >= qubits || target_qubit_begin < 0)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin + qubit_count > qubits)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};
    if (target_qubit_begin <= control_qubit && control_qubit < target_qubit_begin + qubit_count)
        throw std::out_of_range{"Cannot add gate at the given qubit position"};

    //Store circuit array in the host
    af::dim4 dims = internal_circuit.dims();
    std::vector<af::cfloat> vals(dims[0] * dims[1]);
    internal_circuit.host(vals.data());

    // Inserts a given bit to the given value at the given position
    // e.g. insert_bit(5, 1, 1) 5 => b101 = b1011 = 11
    //                                 ^ 
    //                                 1 
    auto insert_bit = [](int value, bool bit, int position) {
        int temp = value & ~ ((1 << position) - 1);
        temp <<= 1;
        temp &= ~(1 << position);
        temp |= (bit ? 1 : 0) << position;
        temp |= value & ((1 << position) - 1);
        return temp;
    };

    //Generates the index of the position of the element in the larget operation circuit matrix from the target qubits of the circuit,
    //the control qubit and the current index being probed of the input circuit matrix
    auto gen_index = [insert_bit](int control, int targets, int target_count, int target_value, int qubit_count, int ii) {
        int index = ii;
        if (control < targets)
        {
            index = insert_bit(index, 1, qubit_count - control - 1 - target_count);
            for (int i = 0; i < target_count; ++i)
                index = insert_bit(index, target_value & (1 << (target_count - 1 - i)), qubit_count - target_count - targets);
        }
        else
        {
            for (int i = 0; i < target_count; ++i)
                index = insert_bit(index, target_value & (1 << (target_count - 1 - i)), qubit_count - target_count - targets - 1);
            index = insert_bit(index, 1, qubit_count - control - 1);
        }
        return index;
    };

    //Create a temporary storage for the operation matrix
    int rem_count = 1 << (qubits - 1 - qubit_count);
    std::vector<af::cfloat> out_temp;
    out_temp.resize(states * states);

    //Generate a identity matrix
    for (int i = 0; i < states; ++i)
        out_temp[i * states + i] = 1.0f;

    //Fill all the entries with the correct values
    for (int i = 0; i < rem_count; ++i)
    {
        for (int m = 0; m < state_count; ++m)
        {
            int ii = gen_index(control_qubit, target_qubit_begin, qubit_count, m, qubits, i);
            for (int n = 0; n < state_count; ++n)
            {
                int jj = gen_index(control_qubit, target_qubit_begin, qubit_count, n, qubits, i);

                //Copy the entries from the input circuit to the correct positions
                out_temp[ii * states + jj] = vals[m * state_count + n];
            }
        }
    }

    //Generate the operation matrix
    af::array gate_mat = af::array(states, states, out_temp.data());

    //Update the circuit matrix
    auto& circuit = qc.circuit();
    circuit = af::matmul(gate_mat, circuit);

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
    temp[0] = state[0] * cosf(angle / 2.f) + state[1] * af::cfloat{ 0.f , -sinf(angle / 2.f) };
    temp[1] = state[0] * af::cfloat{ 0.f , -sinf(angle / 2.f) } + state[1] * cosf(angle / 2.f);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState RotateY_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * cosf(angle / 2.f) + state[1] * -sinf(angle / 2.f);
    temp[1] = state[0] * sinf(angle / 2.f) + state[1] *  cosf(angle / 2.f);

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

QState RotateZ_op(const QState& state, float angle)
{
    af::cfloat temp[2];
    temp[0] = state[0] * af::cfloat{cosf(angle / 2.f) , -sinf(angle / 2.f)} + state[1] * 0.f;
    temp[1] = state[0] * 0.f + state[1] * af::cfloat{cosf(angle / 2.f) , sinf(angle / 2.f)};

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
    temp[1] = state[0] * 0.f  + state[1] * af::cfloat{cosf(angle) , sinf(angle)};

    return { { temp[0].real , temp[0].imag } , { temp[1].real , temp[1].imag} };
}

std::string update_circuit_representation(const std::string& circuit_string, int offset)
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

std::string update_ctrl_circuit_representation(const std::string& circuit_string, int control, int target)
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