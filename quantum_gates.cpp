/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum_gates.h"

namespace aqs
{
   QCircuit Control_GroupPhase(int qubits, int control_qubit, int target_qubit_begin, float angle)
   {
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::invalid_argument{"Invalid control qubit position"};
        if (target_qubit_begin <= control_qubit || target_qubit_begin >= qubits)
            throw std::invalid_argument{"Invalid target qubit begin position"};
    
        QCircuit qc(qubits);
        if (target_qubit_begin == qubits - 1)
        {
            qc << Control_Phase(control_qubit, target_qubit_begin, angle);
        }
        else
        {
            QCircuit cphases(qubits - 1 - target_qubit_begin);
            for (int i = target_qubit_begin; i < qubits; ++i)
                cphases << Phase(i, angle);
        }

        return qc;
   } 

   QCircuit NControl_Gate(int qubits, int control_qubit_begin, int control_qubit_count, int target_qubit_begin, const QCircuit& gate)
   {
        if (gate.qubit_count() >= qubits)
            throw std::invalid_argument{"Gate not supported"};
        if (control_qubit_count <= 0)
            throw std::invalid_argument{"The number of control qubits must be at least 1"};
        if (target_qubit_begin < 0 || (target_qubit_begin + gate.qubit_count()) > qubits)
            throw std::invalid_argument{"Invalid target qubit_begin position"};
        if (control_qubit_begin < 0 || (control_qubit_begin + control_qubit_count) > target_qubit_begin)
            throw std::invalid_argument{"Invalid control_qubit position"};

        QCircuit qc(qubits);
        if (control_qubit_count == 1)
        {
            qc << ControlCircuitGate(gate, control_qubit_begin, target_qubit_begin);
            return qc;
        }
        else
        {
            QCircuit temp(gate.qubit_count() + 1 + target_qubit_begin - control_qubit_begin - control_qubit_count);
            temp << ControlCircuitGate(gate, 0, target_qubit_begin - control_qubit_begin - control_qubit_count + 1);
            for (int i = 0; i < control_qubit_count - 1; ++i)
            {
                QCircuit tmp(temp.qubit_count() + 1);
                tmp << ControlCircuitGate(temp, 0, 1);
                temp = tmp;
            }
            qc << CircuitGate(temp, control_qubit_begin);
            return qc;
        }
   }
}