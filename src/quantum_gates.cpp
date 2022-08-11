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
#include <unordered_set>

namespace aqs
{
    QCircuit Group_Gate(uint32_t qubits, std::vector<uint32_t> target_qubits, QCircuit gate, bool compile)
    {
        if (gate.qubit_count() != 1)
            throw std::invalid_argument{"Gate not supported"};

        std::sort(target_qubits.begin(), target_qubits.end());
        if (target_qubits.size() > qubits)
            throw std::invalid_argument{"Cannot add more target qubits than there are total qubits"};

        gate.compile();
        QCircuit qc{qubits};

        /*
        af::array circuit = af::identity(fast_pow2(target_qubits.front()), fast_pow2(target_qubits.front()), c32);
        circuit = tensor_product(circuit, gate.circuit());

        for (std::size_t i = 1; i < target_qubits.size(); ++i)
        {
            const auto& curr = target_qubits[i];
            const auto& prev = target_qubits[i - 1];
            circuit = tensor_product(circuit,
                      tensor_product(af::identity(fast_pow2(curr - prev - gate.qubit_count()), fast_pow2(curr - prev - gate.qubit_count()), c32), gate.circuit()));
        }

        qc.circuit() = tensor_product(circuit, af::identity(
                                fast_pow2(qubits - 1 - target_qubits.back()), fast_pow2(qubits - 1 - target_qubits.back()), c32));
        */

        for (const auto& qubit : target_qubits)
            qc << Gate{gate, qubit};

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit Control_Group_Gate(uint32_t qubits, uint32_t control_qubit, std::vector<uint32_t> target_qubits,
                                const QCircuit& gate, bool compile)
    {
        if (control_qubit >= qubits)
            throw std::invalid_argument{"Invalid control qubit position"};
        if (gate.qubit_count() != 1)
            throw std::invalid_argument{"Gate not supported"};

        std::sort(target_qubits.begin(), target_qubits.end());
        if (target_qubits.size() >= qubits)
            throw std::invalid_argument{"Cannot add control gate at the given position"};

        auto pivot = std::lower_bound(target_qubits.cbegin(), target_qubits.cend(), control_qubit);
        if (*pivot == control_qubit)
            throw std::invalid_argument{"Cannot add control gate at the target qubit positions"};

        QCircuit qc(qubits);
        auto top_count = std::distance(target_qubits.cbegin(), pivot);
        auto bottom_count = std::distance(pivot, target_qubits.cend());

        if (top_count != 0)
        {
            QCircuit temp(top_count);
            for (auto it = target_qubits.cbegin(); it != pivot; ++it)
                temp << Gate(gate, std::distance(target_qubits.cbegin(), it));

            qc << ControlGate(temp, control_qubit, target_qubits.front());
        }

        if (bottom_count != 0)
        {
            QCircuit temp(bottom_count);
            for (auto it = pivot; it != target_qubits.cend(); ++it)
                temp << Gate(gate, std::distance(pivot, it));

            qc << ControlGate(temp, control_qubit, *pivot);
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit NControl_Gate(uint32_t qubits, uint32_t control_qubit_begin, uint32_t control_qubit_count, uint32_t target_qubit_begin,
                           const QCircuit& gate, bool compile)
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
            qc << ControlGate(gate, control_qubit_begin, target_qubit_begin);
        }
        else
        {
            QCircuit temp(gate.qubit_count() + 1 + target_qubit_begin - control_qubit_begin - control_qubit_count);
            temp << ControlGate(gate, 0, target_qubit_begin - control_qubit_begin - control_qubit_count + 1);
            for (uint32_t i = 0; i < control_qubit_count - 1; ++i)
            {
                QCircuit tmp(temp.qubit_count() + 1);
                tmp << ControlGate(temp, 0, 1);
                temp = std::move(tmp);
            }
            qc << Gate(temp, control_qubit_begin);
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit NControl_Gate(uint32_t qubits, std::vector<uint32_t> control_qubits, uint32_t target_qubit_begin,
                           const QCircuit& gate, bool compile)
    {
        if (control_qubits.size() == 0)
            throw std::invalid_argument{"Number of control qubits must be at least one"};
        if (gate.qubit_count() + control_qubits.size() > qubits)
            throw std::invalid_argument{"Invalid number of qubits for number of control of qubits and gate count"};
        if (target_qubit_begin + gate.qubit_count() > qubits)
            throw std::invalid_argument{"Invalid target qubit begin position"};

        std::sort(control_qubits.begin(), control_qubits.end());
        if (control_qubits.back() >= qubits)
            throw std::invalid_argument{"Cannot add control gate at the given position"};

        auto pivot = std::lower_bound(control_qubits.begin(), control_qubits.end(), target_qubit_begin);
        if (pivot != control_qubits.end() && *pivot < target_qubit_begin + gate.qubit_count())
            throw std::invalid_argument{"Cannot add control gate at the target qubit positions"};

        std::vector<uint32_t> top(control_qubits.begin(), pivot);
        std::vector<uint32_t> bottom(pivot, control_qubits.end());

        QCircuit current = gate;

        // Fit the gate to the fit up to the bottom of the target qubit circuit
        if (bottom.size() != 0)
        {
            // Add bottom control qubits
            for (auto it = bottom.begin(); it != bottom.end() - 1; ++it)
            {
                const auto& control_qubit = *it;
                QCircuit temp(control_qubit - target_qubit_begin + 1);
                temp << ControlGate(current, control_qubit - target_qubit_begin, 0);
                current = std::move(temp);
            }

            QCircuit temp(qubits - target_qubit_begin);
            temp << ControlGate(current, bottom.back() - target_qubit_begin, 0);
            current = std::move(temp);
        }

        // Fit the gate to the fit up to the top of the target qubit circuit
        if (top.size() != 0)
        {
            // Add top control qubits
            uint32_t prev = target_qubit_begin;
            for (auto it = top.rbegin(); it != top.rend() - 1; ++it)
            {
                const auto& control_qubit = *it;
                QCircuit temp(qubits - control_qubit);
                temp << ControlGate{ current , 0 , prev - control_qubit };
                prev = control_qubit;
                current = std::move(temp);
            }

            QCircuit temp(qubits);
            temp << ControlGate(current, top.front(), qubits - current.qubit_count());
            current = std::move(temp);
        }

        if (compile)
            current.compile();

        return current;
   }

   QCircuit Rewire_Gate(uint32_t qubits, const std::vector<uint32_t>& new_qubit_positions, const QCircuit& gate, bool compile)
   {
        if (new_qubit_positions.size() != gate.qubit_count())
            throw std::invalid_argument{"New qubit positions must map all the qubits in the gate"};
        if (gate.qubit_count() > qubits)
            throw std::domain_error{"Cannot rewire circuit to a lower number of qubits"};

        std::vector<uint32_t> sort_positions = new_qubit_positions;
        std::partial_sort_copy(new_qubit_positions.begin(), new_qubit_positions.end(), sort_positions.begin(), sort_positions.end());
        for (uint32_t i = 1; i < sort_positions.size(); ++i)
        {
            if (sort_positions[i - 1] == sort_positions[i])
                throw std::invalid_argument{"Cannot rewire multiple qubits to the same qubit"};
        }

        std::unordered_set<uint32_t> swapped;
        swapped.reserve(new_qubit_positions.size());

        QCircuit qc{ qubits };
        uint32_t swaps = 0;
        for (uint32_t i = 0; i < new_qubit_positions.size(); ++i)
        {
            if (swapped.find(i) != swapped.end())
                continue;
            swapped.insert(i);

            uint32_t current = i;
            while (i != new_qubit_positions[current])
            {
                qc << Swap{ current , new_qubit_positions[current] };
                current = new_qubit_positions[current];
                swapped.insert(current);

                ++swaps;
            }
        }

        qc << Gate{ gate , 0 };

        for (uint32_t i = 0; i < swaps; ++i)
            qc << *dynamic_cast<Swap*>(qc.gate_list()[swaps - 1 - i].get());

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit Adjoint_Gate(const QCircuit& gate)
    {
        QCircuit qc = gate;
        qc.compile();

        auto& list = qc.gate_list();

        // Reverse the order of the added gates
        std::reverse(list.begin(), list.end());

        // Execute the adjoint of each individual gate
        for (auto& g_ptr : list)
        {
            switch (g_ptr->type())
            {
            case Barrier::static_type():
            case X::static_type():
            case Y::static_type():
            case Z::static_type():
            case H::static_type():
            case Swap::static_type():
            case CSwap::static_type():
            case CX::static_type():
            case CY::static_type():
            case CZ::static_type():
            case CCNot::static_type():
            case Or::static_type():
                break;

            case RotX::static_type():
            {
                auto& g = *dynamic_cast<RotX*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case RotY::static_type():
            {
                auto& g = *dynamic_cast<RotY*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case RotZ::static_type():
            {
                auto& g = *dynamic_cast<RotZ*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case Phase::static_type():
            {
                auto& g = *dynamic_cast<Phase*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case CRotX::static_type():
            {
                auto& g = *dynamic_cast<CRotX*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case CRotY::static_type():
            {
                auto& g = *dynamic_cast<CRotY*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case CRotZ::static_type():
            {
                auto& g = *dynamic_cast<CRotZ*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case CPhase::static_type():
            {
                auto& g = *dynamic_cast<CPhase*>(g_ptr.get());
                g.angle = -g.angle;
                break;
            }
            case Gate::static_type():
            {
                auto& g = *dynamic_cast<Gate*>(g_ptr.get());
                g.internal_circuit = std::make_shared<QCircuit>(Adjoint_Gate(*g.internal_circuit));
                break;
            }
            case ControlGate::static_type():
            {
                auto& g = *dynamic_cast<ControlGate*>(g_ptr.get());
                g.internal_circuit = std::make_shared<QCircuit>(Adjoint_Gate(*g.internal_circuit)); 
                break;
            }
            default:
                throw std::runtime_error{ "Unknown unsupported gate cannot be adjoint" };
            }
        }

        // Reverse the order of the circuit string representation
        const auto& rep = qc.representation();
        std::size_t prev = 0;
        auto mark = rep.find(';', prev);
        std::string current;
        while(mark != std::string::npos)
        {
            current = rep.substr(prev, mark - prev + 1) + current;
            prev = mark + 1;
            mark = rep.find(';', prev);
        }
        qc.representation() = current;

        /*
            Equivalent to

            qc.clear_cache();
            qc.compile();
        */
        // af::transposeInPlace(qc.circuit(), true);
        qc.circuit() = af::transpose(qc.circuit(), true);

        return qc;
    }
}