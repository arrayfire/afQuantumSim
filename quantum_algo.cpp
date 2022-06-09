/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum_algo.h"
#include "quantum_gates.h"

namespace aqs
{
    QCircuit grover_oracle(int search_qubits, int marked_state)
    {
        QCircuit qc(search_qubits);

        for (int i = 0; i < search_qubits; ++i)
        {
            if (!(marked_state & (1 << i)))
                qc << X(i);
        }

        //Generate a N-Control Z gate
        QCircuit z(1); z << Z(0);
        qc << CircuitGate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, z), 0);

        for (int i = 0; i < search_qubits; ++i)
        {
            if (!(marked_state & (1 << i)))
                qc << X(i);
        }
        return qc;
    }

    QCircuit grover_search(int search_qubits, const QCircuit& oracle, int iterations, std::string oracle_name)
    {
        if (oracle.qubit_count() < search_qubits)
            throw std::invalid_argument{"Cannot use given oracle for this qubit circuit"};

        QCircuit qc(oracle.qubit_count());
        for (int i = 0; i < search_qubits; ++i)
            qc << Hadamard(i);

        QCircuit z(1); z << Z(0);
        for (int i = 0; i < iterations; ++i)
        {
            qc << Barrier(false);
            qc << CircuitGate(oracle, 0, oracle_name);
            qc << Barrier(false);

            for (int j = 0; j < search_qubits; ++j)
            {
                qc << Hadamard(j);
                qc << X(j);
            }

            //Generate a N-Control Z gate
            qc << CircuitGate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, z), 0);

            for (int j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << Hadamard(j);
            }
        }

        return qc;
    }

    QCircuit grover_iteration(int search_qubits, const QCircuit& oracle, int iterations)
    {
        QCircuit qc(search_qubits);
        QCircuit z(1); z << Z(0);
        for (int i = 0; i < iterations; ++i)
        {
            qc << CircuitGate(oracle, 0);
            for (int j = 0; j < search_qubits; ++j)
            {
                qc << Hadamard(j);
                qc << X(j);
            }
            qc << CircuitGate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, z), 0);
            for (int j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << Hadamard(j);
            }
        }
        return qc;
    }

    QCircuit fourier_transform(int qubits)
    {
        QCircuit qc(qubits);
        for (int i = qubits - 1; i >= 0; --i)
        {
            qc << Hadamard(i);
            for (int j = 0; j < i; ++j)
                qc << Control_Phase(j, i, aqs::pi / (1 << (i - j)));
        }

        return qc;
    }

    QCircuit inverse_fourier_transform(int qubits)
    {
        QCircuit qc(qubits);
        for (int i = 0; i < qubits; ++i)
        {
            for (int j = i - 1; j >= 0; --j)
                qc << Control_Phase(j, i, -aqs::pi / (1 << (i - j)));

            qc << Hadamard(i);
        }

        return qc;
    }
}