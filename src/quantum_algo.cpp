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

#include <nlopt.hpp>

namespace aqs
{
    QCircuit grover_oracle(uint32_t search_qubits, uint32_t marked_state, bool compile)
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

        if (compile)
            qc.compile();        

        return qc;
    }

    QCircuit grover_search(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations, std::string oracle_name, bool compile)
    {
        if (oracle.qubit_count() < search_qubits)
           throw std::invalid_argument{"Cannot use given oracle for this qubit circuit"};

        QCircuit qc(oracle.qubit_count());
        for (uint32_t i = 0; i < search_qubits; ++i)
            qc << H{i};

        for (uint32_t i = 0; i < iterations; ++i)
        {
            qc << Barrier(false);
            qc << Gate(oracle, 0, oracle_name);
            qc << Barrier(false);
            
            for (uint32_t j = 0; j < search_qubits; ++j)
                qc << H{j};

            //qc << Z{search_qubits - 1};
            for (uint32_t j = 0; j < search_qubits; ++j)
                qc << X{j};
            //qc << Z{search_qubits - 1};
            
            //Generate a N-Control Z gate
            qc << Gate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, Z::gate()), 0);

            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << H(j);
            }
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit grover_iteration(uint32_t search_qubits, const QCircuit& oracle, uint32_t iterations, bool compile)
    {
        QCircuit qc(search_qubits);
        for (uint32_t i = 0; i < iterations; ++i)
        {
            qc << Gate(oracle, 0);
            for (uint32_t j = 0; j < search_qubits; ++j)
                qc << H{j};

            for (uint32_t j = 0; j < search_qubits; ++j)
                qc << X{j};

            qc << Gate(NControl_Gate(search_qubits, 0, search_qubits - 1, search_qubits - 1, Z::gate()), 0);

            for (uint32_t j = 0; j < search_qubits; ++j)
            {
                qc << X(j);
                qc << H(j);
            }
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit fourier_transform(uint32_t qubits, bool compile)
    {
        QCircuit qc(qubits);
        for (int32_t i = qubits - 1; i >= 0; --i)
        {
            qc << H(i);
            for (uint32_t j = 0; j < i; ++j)
                qc << CPhase(j, i, aqs::pi / (1 << (i - j)));
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit inverse_fourier_transform(uint32_t qubits, bool compile)
    {
        QCircuit qc(qubits);
        for (uint32_t i = 0; i < qubits; ++i)
        {
            for (int32_t j = i - 1; j >= 0; --j)
                qc << CPhase(j, i, -aqs::pi / (1 << (i - j)));

            qc << H(i);
        }

        if (compile)
            qc.compile();

        return qc;
    }

    std::vector<std::pair<std::string, af::cfloat>> decompose_hamiltonian(const af::array& hamiltonian, uint32_t qubits)
    {
        std::vector<std::pair<std::string, af::cfloat>> coeffs;

        static af::cfloat ivals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {1.f,0.f} };
        static af::cfloat xvals[] = { {0.f,0.f}, {1.f,0.f}, {1.f,0.f}, {0.f,0.f} };
        static af::cfloat yvals[] = { {0.f,0.f}, {0.f,1.f}, {0.f,-1.f}, {0.f,0.f} };
        static af::cfloat zvals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {-1.f,0.f} };

        static af::array I(2, 2, ivals);
        static af::array X(2, 2, xvals);
        static af::array Y(2, 2, yvals);
        static af::array Z(2, 2, zvals);

        auto hs_product = [](const af::array& lhs, const af::array& rhs) {
            return af::sum<af::cfloat>(af::diag(af::matmul(af::transpose(lhs, true), rhs)));
        };

        std::string str;
        str.reserve(qubits);
        for (uint32_t i = 0; i < fast_pow2(qubits * 2); ++i)
        {
            af::array temp = af::identity(1, 1, c32);
            str.clear();

            for (uint32_t j = 0; j < qubits; ++j)
            {
                int val = (i >> (j * 2)) % 4;
                switch (val)
                {
                case 0:
                    temp = tensor_product(temp, I);
                    str.append("i");
                    break;
                case 1:
                    temp = tensor_product(temp, X);
                    str.append("x");
                    break;
                case 2:
                    temp = tensor_product(temp, Y);
                    str.append("y");
                    break;
                case 3:
                    temp = tensor_product(temp, Z);
                    str.append("z");
                    break;
                }
            }

            auto coeff = af::sum<af::cfloat>(af::diag(af::matmul(af::transpose(temp, true), hamiltonian)));
            if ((coeff != af::cfloat{0.f} || coeff != af::cfloat{-0.f}) && (coeff.real * coeff.real + coeff.imag * coeff.imag > 1e-8))
                coeffs.push_back({str, coeff * (1.f / fast_pow2(qubits))});
        }

        return coeffs;
    }

    af::array compose_hamiltonian(const std::vector<std::pair<std::string, af::cfloat>>& description)
    {
        uint32_t qubits = 0;
        if (!description.empty())
            qubits = description[0].first.length();

        af::array out = af::constant(af::cfloat{}, fast_pow2(qubits), fast_pow2(qubits));
        static af::cfloat ivals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {1.f,0.f} };
        static af::cfloat xvals[] = { {0.f,0.f}, {1.f,0.f}, {1.f,0.f}, {0.f,0.f} };
        static af::cfloat yvals[] = { {0.f,0.f}, {0.f,1.f}, {0.f,-1.f}, {0.f,0.f} };
        static af::cfloat zvals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {-1.f,0.f} };

        static af::array I(2, 2, ivals);
        static af::array X(2, 2, xvals);
        static af::array Y(2, 2, yvals);
        static af::array Z(2, 2, zvals);
        for (const auto& pair : description)
        {
            if (pair.first.length() != qubits)
                throw std::invalid_argument{"Invalid Pauli-description"};

            af::array temp = af::identity(1,1,c32);
            const std::string& str = pair.first;
            auto coeff = pair.second;
            for (char c : str)
            {
                switch (c)
                {
                case 'i':
                    temp = tensor_product(temp, I);
                    break;
                case 'x':
                    temp = tensor_product(temp, X);
                    break;
                case 'y':
                    temp = tensor_product(temp, Y);
                    break;
                case 'z':
                    temp = tensor_product(temp, Z);
                    break;
                }
            }
            out += temp * coeff;
        }

        return out;
    }

    void add_evolved_pauli_product(aqs::QCircuit& qc, const std::string& shape, float coeff)
    {
        if (qc.qubit_count() != shape.length())
            throw std::invalid_argument{"Qubit count and shape length must match"};

        if (shape.find_first_of("xyz") == std::string::npos)
        {
            qc << aqs::Barrier{false};
            for (uint32_t i = 0; i < qc.qubit_count(); ++i)
                qc << aqs::RotZ{i, -coeff} << aqs::Phase{i, coeff};
            qc << aqs::Barrier{false};
            return;
        }

        qc << aqs::Barrier{false};
        uint32_t previous = -1;
        for (uint32_t i = 0; i < qc.qubit_count(); ++i)
        {
            const auto& c = shape[i];
            switch(c)
            {
            case 'x':
                qc << aqs::H{i};
                break;
            case 'y':
                qc << aqs::Y{i} << aqs::RotX{i, -aqs::pi / 2.f};
                break;
            case 'z':
                break;
            default:
                continue;
            }

            if (previous != -1)
                qc << aqs::CX{previous, i};

            previous = i;   
        }

        qc << aqs::RotZ{previous, -coeff * 2.f};

        for (uint32_t i = previous; i > 0; --i)
        {
            const auto& c = shape[i - 1];
            if (c == 'i') continue;

            if (i - 1 != previous)
                qc << aqs::CX{i - 1, previous};

            switch (shape[previous])
            {
            case 'x':
                qc << aqs::H{previous};
                break;
            case 'y':
                qc << aqs::Y{previous} << aqs::RotX{previous, -aqs::pi / 2.f};
                break;
            case 'z':
                break;
            }

            previous = i - 1;
        }

        switch (shape[previous])
        {
        case 'x':
            qc << aqs::H{previous};
            break;
        case 'y':
            qc << aqs::Y{previous} << aqs::RotX{previous, -aqs::pi / 2.f};
            break;
        case 'z':
            break;
        }

        qc << aqs::Barrier{false};
    }

    QCircuit hamiltonian_evolution_circuit(const af::array& hamiltonian, uint32_t steps, bool compile)
    {
        uint32_t qubits = fast_log2(hamiltonian.dims()[0]);
        auto hamiltonian_decomposition = decompose_hamiltonian(hamiltonian, qubits);

        aqs::QCircuit qc{qubits};

        for (const auto& pair : hamiltonian_decomposition)
        {
            const std::string& shape = pair.first;
            af::cfloat coeff = pair.second;

            if (coeff.imag != 0.f)
                throw std::invalid_argument{"Must have a real hamiltonian"};

            add_evolved_pauli_product(qc, shape, coeff.real / (float) steps);
            qc << aqs::Barrier{false};
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit linear_entanglement_varstate(uint32_t qubits, uint32_t depth, const std::vector<float>& values, bool compile)
    {
        if (qubits == 0)
            throw std::invalid_argument{"Number of qubits must be greater than zero"};
        if (depth == 0)
            throw std::invalid_argument{"Depth of circuit must be greater than zero"};
        if (qubits * depth * 2 != values.size())
            throw std::invalid_argument{"Parameter count passed is not the expected received"};

        aqs::QCircuit qc{qubits};

        for (uint32_t i = 0; i < qubits; ++i)
            qc << aqs::RotY{i, values[i]};
        for (uint32_t i = 0; i < qubits; ++i)
            qc << aqs::RotZ{i, values[qubits + i]};

        for (uint32_t j = 0; j < depth - 1; ++j)
        {
            for (uint32_t i = 0; i < qubits - 1; ++i)
                qc << aqs::CX{i, i + 1};

            for (uint32_t i = 0; i < qubits; ++i)
                qc << aqs::RotY{i, values[qubits * (j + 2) + i]};
            for (uint32_t i = 0; i < qubits; ++i)
                qc << aqs::RotZ{i, values[qubits * (j + 3) + i]};
        }

        if (compile)
            qc.compile();

        return qc;
    }

    QCircuit full_entanglement_varstate(uint32_t qubits, uint32_t depth, const std::vector<float>& values, bool compile)
    {
        if (qubits == 0)
            throw std::invalid_argument{"Number of qubits must be greater than zero"};
        if (depth == 0)
            throw std::invalid_argument{"Depth of circuit must be greater than zero"};
        if (qubits * depth * 2 != values.size())
            throw std::invalid_argument{"Parameter count passed is not the expected received"};

        aqs::QCircuit qc{qubits};

        for (uint32_t i = 0; i < qubits; ++i)
            qc << aqs::RotY{i, values[i]};
        for (uint32_t i = 0; i < qubits; ++i)
            qc << aqs::RotZ{i, values[qubits + i]};

        for (uint32_t j = 0; j < depth - 1; ++j)
        {
            for (uint32_t i = 0; i < qubits - 1; ++i)
            {
                for (uint32_t k = i + 1; k < qubits; ++k)
                    qc << aqs::CX{i, k};
            }

            for (uint32_t i = 0; i < qubits; ++i)
                qc << aqs::RotY{i, values[qubits * (j + 2) + i]};
            for (uint32_t i = 0; i < qubits; ++i)
                qc << aqs::RotZ{i, values[qubits * (j + 3) + i]};
        }

        if (compile)
            qc.compile();

        return qc;
    }

    std::pair<float, std::vector<float>> variational_quantum_eigensolver(const af::array& matrix,
                                            float range, VQE state_circuit, float tolerance, uint32_t max_evaluations)
    {
        auto dim = matrix.dims();
        if (dim[0] != dim[1] && fast_pow2(fast_log2(dim[0])) != dim[0] && dim[3] != 1 && dim[4] != 1)
            throw std::invalid_argument{"Cannot solve given matrix"};

        uint32_t qubits = fast_log2(dim[0]);
        uint32_t depth = qubits;
        uint32_t param_count = qubits * depth * 2;

        struct Data_t
        {
            std::vector<float>& param_buff;
            aqs::QCircuit& hamiltonian;
            VQE state;
        };

        auto params = std::vector<double>(param_count);
        af::array random_vals = (af::randu(params.size(), f64) * af::Pi * 2.0) - af::Pi;
        random_vals.host(params.data());

        auto param_buff = std::vector<float>(param_count);
        auto hamil_circuit = hamiltonian_evolution_circuit(matrix, range);
        Data_t data{param_buff, hamil_circuit, state_circuit};

        auto cost_function = [](const std::vector<double>& x, std::vector<double>& gradient, void* void_data)
        {
            Data_t& data = *static_cast<Data_t*>(void_data);
            auto& param_buff = data.param_buff;  
            auto& hamiltonian = data.hamiltonian;
            uint32_t qubits = hamiltonian.qubit_count();

            aqs::QCircuit qc{qubits};
            aqs::QSimulator qs{qubits};

            for (std::size_t i = 0; i < param_buff.size(); ++i)
                param_buff[i] = static_cast<float>(x[i]);

            switch (data.state)
            {
            case VQE::LINEAR:
                qc << aqs::Gate{linear_entanglement_varstate(qubits, qubits, param_buff), 0};
                break;
            case VQE::FULL:
                qc << aqs::Gate{full_entanglement_varstate(qubits, qubits, param_buff), 0};
                break;
            }

            qc.compile();
            qs.generate_statevector();
            qs.simulate(qc);

            auto bra_state = af::transpose(qs.statevector(), true);

            qc << aqs::Gate{hamiltonian, 0};
            qc.compile();
            qs.generate_statevector();
            qs.simulate(qc);


            auto ket_state = qs.statevector();
            af::cfloat expectation = af::matmul(bra_state, ket_state)(0).scalar<af::cfloat>();

            double angle = std::acos(expectation.real);
            if (expectation.imag < 0.f)
                angle = -angle;

            return angle;
        };

        nlopt::opt opt(nlopt::LN_COBYLA, params.size());
        opt.set_min_objective(cost_function, &data);
        opt.set_ftol_rel(tolerance);
        opt.set_maxeval(max_evaluations);

        double result = 0.0;
        opt.optimize(params, result);
        result *= range;

        return { result , param_buff };
    }
}