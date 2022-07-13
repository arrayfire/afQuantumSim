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
namespace aqs
{
    // Search algorithms

    /**
     * @brief Generates a circuit containing a certain number of Grover's Iteration
     *        (Oracle + Phase Amplification)
     * 
     * @param search_qubits number of the qubits of the search space
     * @param oracle the oracle circuit that marks the solutions
     * @param iterations number of iterations of Grover's Iteration
     * @return QCircuit 
     */
    QCircuit grover_iteration(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations);

    /**
     * @brief Generates a circuit which executes Grover' search algorithm using the given oracle and iterations
     * 
     * @param search_qubits number of the qubits of the search space
     * @param oracle the oracle circuit that marks the solutions
     * @param iterations number of iterations of Grover's Iteration
     * @param oracle_name name to assign to the oracle gate
     * @return QCircuit 
     */
    QCircuit grover_search(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations, std::string oracle_name = "");

    /**
     * @brief 
     * 
     * @param search_qubits 
     * @param marked_state 
     * @return QCircuit 
     */
    QCircuit grover_oracle(uint32_t search_qubits, uint32_t marked_state);

    /**
     * @brief Generates the Quantum Fourier Algorithm Circuit for the given amount of qubits
     * 
     * @details Most significant digit of the element in the fourier basis is at the top (0-qubit)
     * 
     * @param qubits number of qubits for QFT workspace
     * @return QCircuit 
     */
    QCircuit fourier_transform(uint32_t qubits);

    /**
     * @brief Generates the Inverse Quantum Fourier Algorithm Circuit for the given amount of qubits
     * 
     * @details Receives the most significant digit of the element in the fourier basis at the top (0-qubit)
     * 
     * @param qubits number of qubits for QFT† workspace
     * @return QCircuit 
     */
    QCircuit inverse_fourier_transform(uint32_t qubits);

    /**
     * @brief Returns the Pauli Product decomposition of a hamiltonian(matrix) using the given number of qubits
     * 
     * @details Each pair contains the Pauli product description and the complex coefficient of the product.
     *          The sum of all Pauli product produces the most accurate description of the hamiltonian given the constraints.
     * 
     * @param hamiltonian matrix to decompose
     * @param qubits number of qubits to use
     * 
     * @return std::vector<std::pair<std::string, af::cfloat>> hamiltonian matrix decomposition description
     */
    std::vector<std::pair<std::string, af::cfloat>> decompose_hamiltonian(const af::array& hamiltonian, uint32_t qubits);

    /**
     * @brief Generates a matrix from the given Pauli-product description
     * 
     * @param description Pauli product description of the matrix
     * 
     * @return af::array generated matrix
     */
    af::array compose_hamiltonian(const std::vector<std::pair<std::string, af::cfloat>>& description);

    /**
     * @brief Generates a circuit representation of the time evolution of the given hamiltonian matrix
     * 
     * @param hamiltonian hamiltonian matrix which will be represented
     * @param steps number of steps the time evolution will be segmented in
     * 
     * @return aqs::QCircuit circuit with the time evolved hamiltonian
     */
    QCircuit hamiltonian_evolution_circuit(const af::array& hamiltonian, uint32_t steps);

    /**
     * @brief Creates a circuit linear variational state generator given the number of qubits, the depth of the circuit,
     *        and the parameter values for the rotation angles.
     * 
     * @param qubits number of qubits of the circuit
     * @param depth depth of the circuit
     * @param values parameters of the rotations
     * 
     * @return aqs::QCircuit circuit of the state
     */
    aqs::QCircuit linear_entanglement_varstate(uint32_t qubits, uint32_t depth, const std::vector<float>& values);

    /**
     * @brief Creates a circuit full variational state generator given the number of qubits, the depth of the circuit,
     *        and the parameter values for the rotation angles.
     * 
     * @param qubits number of qubits of the circuit
     * @param depth depth of the circuit
     * @param values parameters of the rotations
     * 
     * @return aqs::QCircuit circuit of the state
     */
    aqs::QCircuit full_entanglement_varstate(uint32_t qubits, uint32_t depth, const std::vector<float>& values);

    enum class VQE
    {
        LINEAR, FULL
    };

    /**
     * @brief Uses the variational quantum eigensolver method (VQE) to find the smallest eigenvalue for
     *        the given hermitianmatrix
     * 
     * @param matrix hermitian matrix
     * @param range the range where the value will be searched [-range,range]
     * 
     * @note Increasing the value of range will help in assuring the smallest eigenvalue, but may decrease the precision
     *       of the answer.
     *       
     * @return value of the smallest eigenvalue found
     */
    std::pair<float, std::vector<float>> variational_quantum_eigensolver(const af::array& matrix, float range, 
                                          VQE state_circuit = VQE::LINEAR, float tolerance = 0, uint32_t max_evaluations = 100);
}