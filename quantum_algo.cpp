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
    QCircuit grover_oracle(uint32_t search_qubits, uint32_t marked_state)
    {
        if (marked_state >= fast_pow2(search_qubits))
            throw std::invalid_argument{"Marked state should be in the range [0, 2^search_qubits)"};

        QCircuit qc(search_qubits);

        for (uint32_t i = 0; i < search_qubits; ++i)
        {
            if (!(marked_state & (1 << i)))
                qc << X(i);
        }

        //Generate a N-Control Z gate
        qc << Gate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, Z::gate()), 0);

        for (int32_t i = 0; i < search_qubits; ++i)
        {
            if (!(marked_state & (1 << i)))
                qc << X(i);
        }
        return qc;
    }

    QCircuit grover_search(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations, std::string oracle_name)
    {
        if (oracle.qubit_count() < search_qubits)
            throw std::invalid_argument{"Cannot use given oracle for this qubit circuit"};

        QCircuit qc(oracle.qubit_count());
        for (uint32_t i = 0; i < search_qubits; ++i)
            qc << H(i);

        for (uint32_t i = 0; i < iterations; ++i)
        {
            qc << Barrier(false);
            qc << Gate(oracle, 0, oracle_name);
            qc << Barrier(false);

            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << H(j);
                qc << X(j);
            }

            //Generate a N-Control Z gate
            qc << Gate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, Z::gate()), 0);

            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << H(j);
            }
        }

        return qc;
    }

    QCircuit grover_iteration(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations)
    {
        QCircuit qc(search_qubits);
        QCircuit z(1); z << Z(0);
        for (uint32_t i = 0; i < iterations; ++i)
        {
            qc << Gate(oracle, 0);
            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << H(j);
                qc << X(j);
            }
            qc << Gate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, z), 0);
            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << H(j);
            }
        }
        return qc;
    }

    QCircuit fourier_transform(uint32_t qubits)
    {
        QCircuit qc(qubits);
        for (int32_t i = qubits - 1; i >= 0; --i)
        {
            qc << H(i);
            for (uint32_t j = 0; j < i; ++j)
                qc << CPhase(j, i, aqs::pi / (1 << (i - j)));
        }

        return qc;
    }

    QCircuit inverse_fourier_transform(uint32_t qubits)
    {
        QCircuit qc(qubits);
        for (uint32_t i = 0; i < qubits; ++i)
        {
            for (int32_t j = i - 1; j >= 0; --j)
                qc << CPhase(j, i, -aqs::pi / (1 << (i - j)));

            qc << H(i);
        }

        return qc;
    }
}