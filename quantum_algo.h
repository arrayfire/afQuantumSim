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
    QCircuit grover_iteration(int search_qubits, const QCircuit& oracle, int iterations);

    /**
     * @brief Generates a circuit which executes Grover' search algorithm using the given oracle and iterations
     * 
     * @param search_qubits number of the qubits of the search space
     * @param oracle the oracle circuit that marks the solutions
     * @param iterations number of iterations of Grover's Iteration
     * @param oracle_name name to assign to the oracle gate
     * @return QCircuit 
     */
    QCircuit grover_search(int search_qubits, const QCircuit& oracle, int iterations, std::string oracle_name = "");

    /**
     * @brief 
     * 
     * @param search_qubits 
     * @param marked_state 
     * @return QCircuit 
     */
    QCircuit grover_oracle(int search_qubits, int marked_state);

    /**
     * @brief Generates the Quantum Fourier Algorithm Circuit for the given amount of qubits
     * 
     * @details Most significant digit of the element in the fourier basis is at the top (0-qubit)
     * 
     * @param qubits number of qubits for QFT workspace
     * @return QCircuit 
     */
    QCircuit fourier_transform(int qubits);

    /**
     * @brief Generates the Inverse Quantum Fourier Algorithm Circuit for the given amount of qubits
     * 
     * @details Receives the most significant digit of the element in the fourier basis at the top (0-qubit)
     * 
     * @param qubits number of qubits for QFT† workspace
     * @return QCircuit 
     */
    QCircuit inverse_fourier_transform(int qubits);
}