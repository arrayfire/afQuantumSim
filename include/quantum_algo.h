/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#pragma once

#include "quantum.h"

/**
 * @brief ArrayFire Quantum Simulator
 *
 */
namespace aqs {
// Search algorithms

/**
 * @brief Generates a circuit containing a certain number of Grover's Iteration
 *        (Oracle + Phase Amplification)
 *
 * @param search_qubits number of the qubits of the search space
 * @param oracle the oracle circuit that marks the solutions
 * @param iterations number of iterations of Grover's Iteration
 * @param compile flag to compile the circuit before returning
 *
 * @return Grover Iteration QCircuit
 */
QCircuit grover_iteration(uint32_t search_qubits, const QCircuit& oracle,
                          uint32_t iterations, bool compile = false);

/**
 * @brief Generates a circuit which executes Grover' search algorithm using the
 * given oracle and iterations
 *
 * @param search_qubits number of the qubits of the search space
 * @param oracle the oracle circuit that marks the solutions
 * @param iterations number of iterations of Grover's Iteration
 * @param oracle_name name to assign to the oracle gate
 * @param compile flag to compile the circuit before returning
 *
 * @return Grover Search QCircuit
 */
QCircuit grover_search(uint32_t search_qubits, const QCircuit& oracle,
                       uint32_t iterations, std::string oracle_name = "",
                       bool compile = false);

/**
 * @brief Generates a oracle circuit that marks the passed state
 *
 * @details Marks the state by doing a pi-phase shift on the state (flipping the
 * |1> state sign)
 *
 * @param search_qubits number of qubits for the search space
 * @param marked_state state to mark with the oracle
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit
 */
QCircuit grover_oracle(uint32_t search_qubits, uint32_t marked_state,
                       bool compile = false);

/**
 * @brief Generates the Quantum Fourier Algorithm Circuit for the given amount
 * of qubits
 *
 * @details Most significant digit of the element in the fourier basis is at the
 * top (0-qubit)
 *
 * @param qubits number of qubits for QFT workspace
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit
 */
QCircuit fourier_transform(uint32_t qubits, bool compile = false);

/**
 * @brief Generates the Inverse Quantum Fourier Algorithm Circuit for the given
 * amount of qubits
 *
 * @details Receives the most significant digit of the element in the fourier
 * basis at the top (0-qubit)
 *
 * @param qubits number of qubits for QFT† workspace
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit
 */
QCircuit inverse_fourier_transform(uint32_t qubits, bool compile = false);

/**
 * @brief Returns the Pauli Product decomposition of a hamiltonian(matrix) using
 * the given number of qubits
 *
 * @details Each pair contains the Pauli product description and the complex
 * coefficient of the product. The sum of all Pauli product produces the most
 * accurate description of the hamiltonian given the constraints.
 *
 * @param hamiltonian matrix to decompose
 * @param qubits number of qubits to use
 *
 * @return std::vector<std::pair<std::string, af::cfloat>> hamiltonian matrix
 * decomposition description
 */
std::vector<std::pair<std::string, af::cfloat>> decompose_hamiltonian(
    const af::array& hamiltonian, uint32_t qubits);

/**
 * @brief Generates a matrix from the given Pauli-product description
 *
 * @param description Pauli product description of the matrix
 *
 * @return af::array generated matrix
 */
af::array compose_hamiltonian(
    const std::vector<std::pair<std::string, af::cfloat>>& description);

/**
 * @brief Generates a circuit representation of the time evolution of the given
 * hamiltonian matrix
 *
 * @param hamiltonian hamiltonian matrix which will be represented
 * @param steps number of steps the time evolution will be segmented in
 * @param compile flag to compile the circuit before returning
 *
 * @return aqs::QCircuit circuit with the time evolved hamiltonian
 */
QCircuit hamiltonian_evolution_circuit(const af::array& hamiltonian,
                                       uint32_t steps, bool compile = false);

/**
 * @brief Creates a circuit linear variational state generator given the number
 * of qubits, the depth of the circuit, and the parameter values for the
 * rotation angles.
 *
 * @param qubits number of qubits of the circuit
 * @param depth depth of the circuit
 * @param values parameters of the rotations
 * @param compile flag to compile the circuit before returning
 *
 * @return aqs::QCircuit circuit of the state
 */
QCircuit linear_entanglement_varstate(uint32_t qubits, uint32_t depth,
                                      const std::vector<float>& values,
                                      bool compile = false);

/**
 * @brief Creates a circuit full variational state generator given the number of
 * qubits, the depth of the circuit, and the parameter values for the rotation
 * angles.
 *
 * @param qubits number of qubits of the circuit
 * @param depth depth of the circuit
 * @param values parameters of the rotations
 * @param compile flag to compile the circuit before returning
 *
 * @return aqs::QCircuit circuit of the state
 */
QCircuit full_entanglement_varstate(uint32_t qubits, uint32_t depth,
                                    const std::vector<float>& values,
                                    bool compile = false);

enum class VQE { LINEAR, FULL };

/**
 * @brief Uses the variational quantum eigensolver method (VQE) to find the
 * smallest eigenvalue for the given hermitianmatrix
 *
 * @param matrix hermitian matrix
 * @param range the range where the value will be searched [-range,range]
 * @param state_circuit the type of the variational state generator circuit used
 * for the parameters
 * @param tolerance the precision to which the minimum eigenvalue will be search
 * to
 * @param max_evaluations the maximum number of optimizations trials that the
 * algorithm will execute
 *
 * @note Increasing the value of range will help in assuring the smallest
 * eigenvalue, but may decrease the precision of the answer.
 *
 * @return value of the smallest eigenvalue found
 */
std::pair<float, std::vector<float>> variational_quantum_eigensolver(
    const af::array& matrix, float range, VQE state_circuit = VQE::LINEAR,
    float tolerance = 0, uint32_t max_evaluations = 100);
}  // namespace aqs