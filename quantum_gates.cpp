/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum_gates.h"

#include <algorithm>

namespace aqs
{
   QCircuit Control_GroupPhase(uint32_t qubits, uint32_t control_qubit, uint32_t target_qubit_begin, float angle)
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
            for (uint32_t i = target_qubit_begin; i < qubits; ++i)
                cphases << Phase(i, angle);
        }

        return qc;
   } 

   QCircuit NControl_Gate(uint32_t qubits, uint32_t control_qubit_begin, uint32_t control_qubit_count, uint32_t target_qubit_begin, const QCircuit& gate)
   {
        if (gate.qubit_count() >= qubits)
            throw std::invalid_argument{"Gate not supported"};
        if (control_qubit_count == 0)
            throw std::invalid_argument{"The number of control qubits must be at least 1"};
        if ((target_qubit_begin + gate.qubit_count()) > qubits)
            throw std::invalid_argument{"Invalid target qubit_begin position"};
        if ((control_qubit_begin + control_qubit_count) > target_qubit_begin)
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
            for (uint32_t i = 0; i < control_qubit_count - 1; ++i)
            {
                QCircuit tmp(temp.qubit_count() + 1);
                tmp << ControlCircuitGate(temp, 0, 1);
                temp = tmp;
            }
            qc << CircuitGate(temp, control_qubit_begin);
            return qc;
        }
   }

   QCircuit NControl_Gate(uint32_t qubits, std::vector<uint32_t> control_qubits, uint32_t target_qubit_begin, const QCircuit& gate)
   {
        if (control_qubits.size() == 0)
            throw std::invalid_argument{"Number of control qubits must be at least one"};
        if (gate.qubit_count() + control_qubits.size() > qubits)
            throw std::invalid_argument{"Invalid number of qubits for number of control of qubits and gate count"};
        if (target_qubit_begin + gate.qubit_count() > qubits)
            throw std::invalid_argument{"Invalid target qubit begin position"};

        std::sort(control_qubits.begin(), control_qubits.end());
        if (control_qubits.back() <= target_qubit_begin + gate.qubit_count())
            throw std::invalid_argument{"Cannot add control gate at the target qubit positions"};

        auto pivot = std::lower_bound(control_qubits.begin(), control_qubits.end(), target_qubit_begin);
        if (pivot != control_qubits.end() && *pivot < target_qubit_begin + gate.qubit_count())
            throw std::invalid_argument{"Cannot add control gate at the target qubit positions"};

        std::vector<uint32_t> top(control_qubits.begin(), pivot);
        std::vector<uint32_t> bottom(pivot, control_qubits.end());

        QCircuit qc(qubits);
        QCircuit current = gate;

        // Add bottom control qubits
        for (auto it = bottom.begin(); it < bottom.end(); ++it)
        {
            const auto& control_qubit = *it;
            QCircuit temp(control_qubit - target_qubit_begin + 1);
            temp << ControlCircuitGate(current, control_qubit, 0);
            current = std::move(temp);
        }
        // Fit the gate to the fit up to the bottom of the target qubit circuit
        if (bottom.size() != 0)
        {
            QCircuit temp(qubits - target_qubit_begin);
            temp << ControlCircuitGate(current, bottom.back(), 0);
            current = std::move(temp);
        }

        // Add top control qubits
        for (auto it = top.rbegin(); it < top.rend(); ++it)
        {
            const auto& control_qubit = *it;
            QCircuit temp(target_qubit_begin - control_qubit + current.qubit_count());
            temp << ControlCircuitGate(current, 0, target_qubit_begin - control_qubit);
            current = std::move(temp);
        }
        // Fit the gate to the fit up to the top of the target qubit circuit
        if (top.size() != 0)
        {
            QCircuit temp(qubits);
            temp << ControlCircuitGate(current, top.front(), qubits - current.qubit_count());
            current = std::move(temp);
        }

        return qc;
   }
}