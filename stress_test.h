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

        //Find the n-qubit X matrix for the target qubit as I_L @ Y @ I_R where @ is the tensor product
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

        //Find the n-qubit X matrix for the target qubit as I_L @ Z @ I_R where @ is the tensor product
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

        //Find the n-qubit X matrix for the target qubit as I_L @ Phase @ I_R where @ is the tensor product
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

class Hadamard_old : public aqs::QGate
{
public:
    Hadamard_old(uint32_t qubit_) : target_qubit{qubit_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        af::array left_identity = af::identity(fast_pow2(target_qubit), fast_pow2(target_qubit), c32);
        af::array right_identity = af::identity(fast_pow2(qubits - target_qubit - 1), fast_pow2(qubits - target_qubit - 1), c32);

        //Find the n-qubit X matrix for the target qubit as I_L @ H @ I_R where @ is the tensor product
        af::array temp = tensor_product(left_identity, hadamard_matrix());
        temp = tensor_product(temp, right_identity);

        auto& circuit = qc.circuit();
        circuit = af::matmul(temp, circuit);
        return qc;
    }
    std::string to_string() const override { return "H," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};

    static const af::array& hadamard_matrix()
    {
        static af::cfloat vals[] = {
            {1.f, 0.f} , {0.f,1.f},
            {0.f, 0.f} , {-1.f, 0.f}
        };
        static af::array mat{2, 2, vals};
        return mat;
    }
};

class Control_X_old : public aqs::QGate
{
public:
    Control_X_old(uint32_t control, uint32_t target) : control_qubit{control} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit < 0 || target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);

        int control_mask = 1 << (qubits - 1 - control_qubit);
        int target_mask = 1 << (qubits - 1 - target_qubit);

        //Generate the sparse matrix positions
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            rows[i + 1] = i + 1;
            cols[i] = (i & control_mask) ? (i ^ target_mask) : i;
        }

        //All entries are (1.f, 0.0f), generate the sparse matrix
        af::array cx_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f}, states), af::array(states + 1, rows.data()), af::array(states, cols.data()));

        //Update the circuit matrix by matrix multiplication
        auto& circuit = qc.circuit();
        circuit = af::matmul(cx_matrix, circuit);

        return qc;
    }
    std::string to_string() const override { return "CX,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};
    uint32_t control_qubit{};
};

class Control_Y_old : public aqs::QGate
{
public:
    Control_Y_old(uint32_t control, uint32_t target) : control_qubit{control} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit < 0 || target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);
        std::vector<af::cfloat> vals(states);

        int control_mask = 1 << (qubits - 1 - control_qubit);
        int target_mask = 1 << (qubits - 1 - target_qubit);

        // Generate the sparse matrix entry positions and values
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            vals[i] = (i & control_mask) ? (i & target_mask ? af::cfloat{0.f, 1.f} : af::cfloat{0.f, -1.f}) : 1.0f;
            rows[i + 1] = i + 1;
            cols[i] = (i & control_mask) ? (i ^ target_mask) : i;
        }

        //Generate operation matrix
        af::array cy_matrix = af::sparse(states, states, af::array(states, vals.data()), af::array(states + 1, rows.data()),
                                        af::array(states, cols.data()));

        //Update the circuit matrix by matrix multiplication
        auto& circuit = qc.circuit();
        circuit = af::matmul(cy_matrix, circuit);

        return qc;
    }
    std::string to_string() const override { return "CY,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};
    uint32_t control_qubit{};
};

class Control_Z_old : public aqs::QGate
{
public:
    Control_Z_old(uint32_t control, uint32_t target) : control_qubit{control} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit < 0 || target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);
        std::vector<af::cfloat> vals(states);

        int control_mask = 1 << (qubits - 1 - control_qubit);
        int target_mask = 1 << (qubits - 1 - target_qubit);
        int mask = control_mask | target_mask;

        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            vals[i] = (i & mask) != mask ? 1.0f : -1.0f;
            rows[i + 1] = i + 1;
            cols[i] = i;
        }

        af::array cz_matrix = af::sparse(states, states, af::array(states, vals.data()),
                                        af::array(states + 1, rows.data()), af::array(states, cols.data()));

        auto& circuit = qc.circuit();
        circuit = af::matmul(cz_matrix, circuit);

        return qc;
    }
    std::string to_string() const override { return "CZ,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};
    uint32_t control_qubit{};
};

class Control_Phase_old : public aqs::QGate
{
public:
    Control_Phase_old(uint32_t control, uint32_t target, float angle_) : control_qubit{control} , target_qubit{target} , angle{angle_} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit < 0 || target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);
        std::vector<af::cfloat> vals(states);

        int control_mask = 1 << (qubits - 1 - control_qubit);
        int target_mask = 1 << (qubits - 1 - target_qubit);
        int mask = control_mask | target_mask;

        //Generate the sparse entries locations and values
        af::cfloat rot = {cosf(angle), sinf(angle)};
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            vals[i] = ((i & mask) == mask)? rot : 1.0f;
            rows[i + 1] = i + 1;
            cols[i] = i;
        }

        //Generate the operation matrix
        af::array cphase_matrix = af::sparse(states, states, af::array(states, vals.data()),
                                        af::array(states + 1, rows.data()), af::array(states, cols.data()));

        //Update the circuit matrix by matrix multiplication
        auto& circuit = qc.circuit();
        circuit = af::matmul(cphase_matrix, circuit);

        return qc;
    }
    std::string to_string() const override { return "Phase," + std::to_string(target_qubit) + ";"; }
    uint32_t control_qubit{};
    uint32_t target_qubit{};
    float angle{};
};

class CControl_Not_old : public aqs::QGate
{
public:
    CControl_Not_old(uint32_t controlA, uint32_t controlB, uint32_t target) : control_qubit_A{controlA} ,
        control_qubit_B{controlB} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();

        if (qubits < 3 || control_qubit_A == target_qubit || control_qubit_B == target_qubit)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit_A >= qubits || control_qubit_A < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (control_qubit_B >= qubits || control_qubit_B < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit >= qubits || target_qubit < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);

        int control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
        int target_mask = 1 << (qubits - 1 - target_qubit);

        //Generate the matrix entry indices
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            rows[i + 1] = i + 1;
            cols[i] = (i & control_mask) == control_mask ? (i ^ target_mask) : i;
        }

        //All entries are 1, generate the operation matrix 
        af::array ccnot_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                            af::array(states + 1, rows.data()), af::array(states, cols.data()));

        //Update the circuit matrix by matrix multiplication
        auto& circuit = qc.circuit();
        circuit = af::matmul(ccnot_matrix, circuit);

        return qc;
    }

    std::string to_string() const override { return "CCNot,2,1:" + std::to_string(control_qubit_A) + "," + std::to_string(control_qubit_B) +
                                             "," + std::to_string(target_qubit) + ";"; }
    uint32_t control_qubit_A{};
    uint32_t control_qubit_B{};
    uint32_t target_qubit{};
};

class Or_old : public aqs::QGate
{
public:
    Or_old(uint32_t controlA, uint32_t controlB, uint32_t target) : control_qubit_A{controlA} ,
        control_qubit_B{controlB} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (control_qubit_A >= qubits || control_qubit_A < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (control_qubit_B >= qubits || control_qubit_B < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit >= qubits || target_qubit < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);

        int control_mask = (1 << (qubits - 1 - control_qubit_A)) | (1 << (qubits - 1 - control_qubit_B));
        int target_mask = 1 << (qubits - 1 - target_qubit);

        //Generate the sparse matrix entries
        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            rows[i + 1] = i + 1;
            cols[i] = (i & control_mask) ? i ^ target_mask : i;
        }

        //Generate the operation matrix
        af::array or_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                        af::array(states + 1, rows.data()), af::array(states, cols.data()));

        //Update the circuit matrix by matrix multiplication
        auto& circuit = qc.circuit();
        circuit = af::matmul(or_matrix, circuit);

        return qc;
    }

    std::string to_string() const override { return "Or,2,1:" + std::to_string(control_qubit_A) + "," + std::to_string(control_qubit_B) +
                                             "," + std::to_string(target_qubit) + ";"; }
    uint32_t control_qubit_A{};
    uint32_t control_qubit_B{};
    uint32_t target_qubit{};
};

class Control_Swap_old : public aqs::QGate
{
public:
    Control_Swap_old(uint32_t control, uint32_t targetA, uint32_t targetB) : control_qubit{control} , target_qubit_A{targetA} ,
        target_qubit_B{targetB} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        if (qubits < 3)
            throw std::domain_error{"Gate not supported for given circuit"};
        if (target_qubit_A >= qubits || target_qubit_A < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit_B >= qubits || target_qubit_B < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (control_qubit >= qubits || control_qubit < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        std::vector<int> cols(states), rows(states + 1);

        int control_mask = 1 << (qubits - 1 - control_qubit);
        int target_maskA = 1 << (qubits - 1 - target_qubit_A);
        int target_maskB = 1 << (qubits - 1 - target_qubit_B);
        int target_mask = target_maskA | target_maskB;

        rows[0] = 0;
        for (int i = 0; i < states; ++i)
        {
            rows[i + 1] = i + 1;
            if ((i & control_mask) && ((i & target_mask) == target_maskA || (i & target_mask) == target_maskB))
            {
                cols[i] = (i ^ target_maskA) ^ target_maskB;
            }
            else
                cols[i] = i;
        }

        af::array cswap_matrix = af::sparse(states, states, af::constant(af::cfloat{ 1.0f , 0.0f }, states),
                                af::array(states + 1, rows.data()), af::array(states, cols.data()));

        auto& circuit = qc.circuit();
        circuit = af::matmul(cswap_matrix, circuit);

        return qc;
    }

    std::string to_string() const override { return "CSwap,1,2:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit_A) +
                                             "," + std::to_string(target_qubit_B) + ";"; }
    uint32_t control_qubit{};
    uint32_t target_qubit_A{};
    uint32_t target_qubit_B{};
};

class Control_Hadamard_old : public aqs::QGate
{
public:
    Control_Hadamard_old(uint32_t control, uint32_t target) : control_qubit{control} , target_qubit{target} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        if (qubits < 2)
            throw std::domain_error{"Gate not supported for given simulation"};
        if (control_qubit < 0 || control_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit < 0 || target_qubit >= qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        aqs::QCircuit h(1);
        h << aqs::Hadamard{0};

        qc << aqs::ControlCircuitGate{h, control_qubit, target_qubit};

        return qc;
    }
    std::string to_string() const override { return "CHadamard,1,1:" + std::to_string(control_qubit) + "," + std::to_string(target_qubit) + ";"; }
    uint32_t target_qubit{};
    uint32_t control_qubit{};
};

class CircuitGate_old : public aqs::QGate
{
public:
    CircuitGate_old(const aqs::QCircuit& circuit_, uint32_t target_qubit_begin_, std::string name = "")
        : internal_circuit{circuit_.circuit()}, qubit_count{circuit_.qubit_count()}, target_qubit_begin{target_qubit_begin_},
            representation{name} {}
    aqs::QCircuit& operator()(aqs::QCircuit& qc) const override
    {
        const int qubits = qc.qubit_count();
        const int states = qc.state_count();
        const int state_count = fast_pow2(qubit_count);
        if (target_qubit_begin >= qubits || target_qubit_begin < 0)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};
        if (target_qubit_begin + qubit_count > qubits)
            throw std::out_of_range{"Cannot add gate at the given qubit position"};

        //Store circuit array in the host
        af::dim4 dims = internal_circuit.dims();
        std::vector<af::cfloat> vals(dims[0] * dims[1]);
        internal_circuit.host(vals.data());

        // Inserts a given bit to the given value at the given position
        // e.g. insert_bit(5, 1, 1) 5 => b101 = b1011 = 11
        //                                 ^ 
        //                                 1 
        auto insert_bit = [](int value, bool bit, int position) {
            int temp = value & ~ ((1 << position) - 1);
            temp <<= 1;
            temp &= ~(1 << position);
            temp |= (bit ? 1 : 0) << position;
            temp |= value & ((1 << position) - 1);
            return temp;
        };

        //Generates the index of the position of the element in the larget operation circuit matrix from the target qubits of the circuit,
        //and the current index being probed of the input circuit matrix
        auto gen_index = [insert_bit](int targets, int target_count, int target_value, int qubit_count, int ii) {
            int index = ii;
            for (int i = 0; i < target_count; ++i)
                index = insert_bit(index, target_value & (1 << (target_count - 1 - i)), qubit_count - target_count - targets);
            return index;
        };

        //Create a temporary storage for the operation matrix
        int rem_count = fast_pow2(qubits - qubit_count);
        std::vector<af::cfloat> out_temp;
        out_temp.resize(states * states);

        //Generate a identity matrix
        for (int i = 0; i < states; ++i)
            out_temp[i * states + i] = 1.0f;

        //Fill all the entries with the correct values
        for (int i = 0; i < rem_count; ++i)
        {
            for (int m = 0; m < state_count; ++m)
            {
                int ii = gen_index(target_qubit_begin, qubit_count, m, qubits, i);
                for (int n = 0; n < state_count; ++n)
                {
                    int jj = gen_index(target_qubit_begin, qubit_count, n, qubits, i);

                    //Copy the entries from the input circuit to the correct positions
                    out_temp[ii * states + jj] = vals[m * state_count + n];
                }
            }
        }

        //Generate the operation matrix
        af::array gate_mat = af::array(states, states, out_temp.data());

        //Update the circuit matrix
        auto& circuit = qc.circuit();
        circuit = af::matmul(gate_mat, circuit);

        return qc;
    }
    std::string to_string() const override { return representation; }
    af::array internal_circuit;
    std::string representation;
    uint32_t qubit_count;
    uint32_t target_qubit_begin;
};