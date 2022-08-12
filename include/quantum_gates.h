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

namespace aqs {
/**
 * @brief Creates circuit containing the same gate added at multiple target
 * qubits
 *
 * @param qubits number of qubits in the circuit
 * @param target_qubits list of the qubits where the gate will be added
 * @param gate gate to be added
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit circuit containing the group of gates
 */
QCircuit Group_Gate(uint32_t qubits, std::vector<uint32_t> target_qubits,
                    QCircuit gate, bool compile = false);

/**
 * @brief Creates a circuit containing a group of gates controlled by a qubit
 *
 * @param qubits number of qubits in the circuit
 * @param control_qubit position of the qubit controlling the gates
 * @param target_qubits list of the qubits where the gate will be added
 * @param gate gate to be added
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit circuit containing the controlled gates
 */
QCircuit Control_Group_Gate(uint32_t qubits, uint32_t control_qubit,
                            std::vector<uint32_t> target_qubits,
                            const QCircuit& gate, bool compile = false);

/**
 * @brief Creates a circuit containing a gate controlled by a given number of
 * consecutive qubits
 *
 * @param qubits number of qubits in the circuit
 * @param control_qubit_begin starting position where the qubits will be added
 * @param control_qubit_count number of control qubits to be added consecutively
 * @param target_qubit_begin position to add the gate to
 * @param gate gate to be added
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit circuit with the gate controlled by N consecutive qubits
 */
QCircuit NControl_Gate(uint32_t qubits, uint32_t control_qubit_begin,
                       uint32_t control_qubit_count,
                       uint32_t target_qubit_begin, const QCircuit& gate,
                       bool compile = false);

/**
 * @brief Creates a circuit containing a gate controlled by the given number of
 * qubits
 *
 * @param qubits number of qubits in the circuit
 * @param control_qubits list of the control qubits which will be set
 * @param target_qubit_begin the position where the gate will be added
 * @param gate gate to be added
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit circuit with the gate controlled by N qubits
 */
QCircuit NControl_Gate(uint32_t qubits, std::vector<uint32_t> control_qubits,
                       uint32_t target_qubit_begin, const QCircuit& gate,
                       bool compile = false);

/**
 * @brief Creates a circuit where the connections have been rewired for the
 * qubits
 *
 * @note the new qubit positions must map every single qubit in the gate to a
 * different qubit desired in the circuit
 *
 * @param qubits number of qubits in the circuit
 * @param new_qubit_positions the new positions that the qubits will be mapped
 * @param gate gate to be added
 * @param compile flag to compile the circuit before returning
 *
 * @return QCircuit rewired circuit
 */
QCircuit Rewire_Gate(uint32_t qubits,
                     const std::vector<uint32_t>& new_qubit_positions,
                     const QCircuit& gate, bool compile = false);

/**
 * @brief Returns the adjoint (inverse/reversed) version of the gate
 *
 * @param gate gate to be adjoint
 * @return QCircuit adjoint gate
 */
QCircuit Adjoint_Gate(const QCircuit& gate);
}  // namespace aqs