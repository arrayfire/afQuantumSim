#pragma once

#include "quantum.h"

class X_old : public aqs::QGate
{
public:
    X_old(uint32_t qubit_) : target_qubit{qubit_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
        af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

        //Find the n-qubit X matrix for the target qubit as I_L @ X @ I_R where @ is the tensor product
        af::array temp = tensor_product(left_identity, x_matrix());
        temp = tensor_product(temp, right_identity);

        auto& circuit = qc.circuit();
        circuit = af::matmul(temp, circuit);
        return qc;
    }
    std::string to_string() const override { return "X," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};

    static const af::array& x_matrix()
    {
        static af::cfloat vals[] = {
            {0.f, 0.f} , {1.f, 0.f},
            {1.f, 0.f} , {0.f, 0.f}
        };
        static af::array mat{2, 2, vals};
        return mat;
    }
};

class Y_old : public aqs::QGate
{
public:
    Y_old(uint32_t qubit_) : target_qubit{qubit_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
        af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

        //Find the n-qubit X matrix for the target qubit as I_L @ X @ I_R where @ is the tensor product
        af::array temp = tensor_product(left_identity, y_matrix());
        temp = tensor_product(temp, right_identity);

        auto& circuit = qc.circuit();
        circuit = af::matmul(temp, circuit);
        return qc;
    }
    std::string to_string() const override { return "Y," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};

    static const af::array& y_matrix()
    {
        static af::cfloat vals[] = {
            {0.f, 0.f} , {0.f,-1.f},
            {0.f, 1.f} , {0.f, 0.f}
        };
        static af::array mat{2, 2, vals};
        return mat;
    }
};

class Z_old : public aqs::QGate
{
public:
    Z_old(uint32_t qubit_) : target_qubit{qubit_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
        af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

        //Find the n-qubit X matrix for the target qubit as I_L @ X @ I_R where @ is the tensor product
        af::array temp = tensor_product(left_identity, z_matrix());
        temp = tensor_product(temp, right_identity);

        auto& circuit = qc.circuit();
        circuit = af::matmul(temp, circuit);
        return qc;
    }
    std::string to_string() const override { return "Z," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};

    static const af::array& z_matrix()
    {
        static af::cfloat vals[] = {
            {1.f, 0.f} , {0.f,1.f},
            {0.f, 0.f} , {-1.f, 0.f}
        };
        static af::array mat{2, 2, vals};
        return mat;
    }
};

class Phase_old : public aqs::QGate
{
public:
    Phase_old(uint32_t qubit_, float angle_) : target_qubit{qubit_} , angle{angle_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
        af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

        //Find the n-qubit X matrix for the target qubit as I_L @ X @ I_R where @ is the tensor product
        const af::cfloat vals[] = {
            { 1.0f, 0.f } , { 0.f , 0.f },
            { 0.f , 0.f } , { std::cos(angle) , std::sin(angle) }
        };
        af::array phase_matrix(2, 2, vals);
        af::array temp = tensor_product(left_identity, phase_matrix);
        temp = tensor_product(temp, right_identity);

        auto& circuit = qc.circuit();
        circuit = af::matmul(temp, circuit);
        return qc;
    }
    std::string to_string() const override { return "Phase," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};
    float angle{};
};

class Swap_old : public aqs::QGate
{
public:
    Swap_old(int target_qubit_A_, int target_qubit_B_) : target_qubit_A{target_qubit_A_} , target_qubit_B{target_qubit_B_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
            const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (target_qubit_A < 0 || target_qubit_A >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit_B < 0 || target_qubit_B >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);

        int control_mask = 1 << (qubits - 1 - target_qubit_A);
        int target_mask = 1 << (qubits - 1 - target_qubit_B);
        int mask = control_mask | target_mask;

        //Generate all entry indices
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            rows[i + 1] = i + 1;
            if ((i & mask) == control_mask || (i & mask) == target_mask)
                cols[i] = (i ^ target_mask) ^ control_mask;
            else
                cols[i] = i;
        }

        //All entries are (1.0f, 0.0f), generate operation matrix from the generated indices
        af::array swap_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                        af::array(states + 1, rows.data()), af::array(states, cols.data()));
        
        auto& circuit = qc.circuit();
        circuit = af::matmul(swap_matrix, circuit);

        return qc;
    }

    std::string to_string() const override
    {
        return "Swap," + std::to_string(target_qubit_A) + "," + std::to_string(target_qubit_B) + ";";
    }
    
    int target_qubit_A{}, target_qubit_B{};
};

