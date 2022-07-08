#include <iostream>
#include <vector>
#include <numeric>
#include <random>

#include <nlopt.hpp>
#include "quantum.h"
#include "quantum_visuals.h"

template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vec)
{
    stream << "[ ";
    for (const auto& val : vec)
        stream << val << " ";
    stream << "]\n";
    return stream;
}

template<typename T, typename U>
std::ostream& operator<<(std::ostream& stream, const std::vector<std::pair<T, U>>& vec)
{
    stream << "[ ";
    for (const auto& val : vec)
        stream << val.first << " ";
    stream << "]\n" << "[ ";
    for (const auto& val : vec)
        stream << val.second << " ";
    stream << "]\n";
    return stream;
}

std::vector<double> matrix_vec_mult(const std::vector<double>& mat, const std::vector<double>& vec)
{
    if (mat.size() % vec.size() != 0)
        throw std::invalid_argument{"Invalid multiplication"};

    std::vector<double> out(mat.size() / vec.size(), 0.0);
    for (std::size_t i = 0; i < out.size(); ++i)
    {
        for (std::size_t j = 0; j < vec.size(); ++j)
            out[i] += mat[i * out.size() + j] * vec[j];
    }

    return out;
}

double abs_error(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    std::vector<double> mat = {
        1.0, 1.0,
        1.0,-1.0
    };

    std::vector<double> vec = {
        7.f, 3.f
    };

    auto res = matrix_vec_mult(mat, x);

    /*return std::accumulate(x.cbegin(), x.cend(), 0.0, [](double a, double b) {
        return std::abs(a) + std::abs(b);
    });*/
    double error = 0;
    for (std::size_t i = 0; i < res.size(); ++i)
    {
        double diff = res[i] - vec[i];
        error += diff * diff;
    }

    return error;
}

double constraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    return std::accumulate(x.cbegin(), x.cend(), 0.0);
}

static std::vector<double> target;
static std::unique_ptr<aqs::QCircuit> qc_ptr;
static std::unique_ptr<aqs::QSimulator> qs_ptr;
double obj_cost(const std::vector<double>& params, std::vector<double>& grad, void* data)
{
    aqs::QCircuit& qc = *qc_ptr;
    aqs::QSimulator& qs = *qs_ptr;

    const auto& lamda = params[0];
    const auto& theta = params[1];
    const auto& phi = params[2];

    qc.reset_circuit();
    qc << aqs::Phase{0, (float)lamda} << aqs::RotY{0, (float)theta} << aqs::Phase{0, (float)phi};
    qc.generate_circuit();

    qs.generate_global_state();
    qs.simulate(qc);

    /*
    auto measures = qs.probabilities();

    double cost = 0.0;
    for (int i = 0; i < 3; ++i)
        cost += std::abs(measures[i] - target[i]);

    return cost;
    */
    af::array bra_state = af::transpose(qs.global_state(), true);
    aqs::QCircuit q(1); q << aqs::Phase{0, (float)lamda} << aqs::RotY{0, (float)theta} << aqs::Phase{0, (float)phi} << aqs::Z{0};
    q.generate_circuit();
    aqs::QSimulator s(1); s.simulate(q);
    af::array ket_state = s.global_state();
    double expect = static_cast<double>(af::real(af::matmul(bra_state, ket_state))(0).scalar<float>());

    return expect;
}

static std::random_device dv{};
static std::mt19937 rd{dv()};
static std::uniform_real_distribution<double> dist{0.0, 1.0};
void begin()
{
    qc_ptr = std::make_unique<aqs::QCircuit>(1);
    qs_ptr = std::make_unique<aqs::QSimulator>(1);

    target = std::vector<double>(2);
    target[0] = dist(rd);
    target[1] = 1.0 - target[0];

    std::cout << target << std::endl;
}

void optimize()
{
    nlopt::opt opt{nlopt::LN_COBYLA, 3};

    opt.set_min_objective((nlopt::vfunc)obj_cost, nullptr);

    opt.set_xtol_rel(1e-5);
    opt.set_ftol_rel(1e-5);
    opt.set_maxeval(10000);
    opt.set_maxtime(20);
    opt.set_initial_step(1e-4);

    double pi = af::Pi;
    opt.set_lower_bounds({-pi, -pi, -pi});
    opt.set_upper_bounds({pi, pi, pi});

    std::uniform_real_distribution<double> d{-pi, pi};
    std::vector<double> out = { d(rd) , d(rd), d(rd) };
    double val = 0.0;
    opt.optimize(out, val);

    std::cout << val << std::endl;
    std::cout << out;
    std::cout << target;
    std::cout << qs_ptr->probabilities();
}

aqs::QCircuit linear_vqe(uint32_t qubits, uint32_t depth, const std::vector<float>& values)
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

    qc.generate_circuit();
    return qc;
}

aqs::QCircuit hamiltonian(uint32_t qubits)
{
    aqs::QCircuit qc{qubits};

    for (uint32_t i = 0; i < qubits; ++i)
        qc << aqs::H{i};

    qc.generate_circuit();
    return qc;
}

struct Data
{
    uint32_t qubits;
};

std::vector<float> get_coefficients(uint32_t qubits, const af::array& hamiltonian)
{
    std::vector<float> coeff;
    coeff.reserve(qubits * 4);
    aqs::QCircuit qc(qubits);

    auto hilsch_product = [](const af::array& lhs, const af::array& rhs) {
        return af::sum<float>(af::diag(af::real(af::matmul(af::transpose(lhs, true), rhs))));
    };

    for (uint32_t i = 0; i < qubits; ++i)
    {
        af::array identity = af::identity(qc.state_count(), qc.state_count(), c32);
        qc << aqs::X{i};
        qc.generate_circuit();
        af::array pauli_x = qc.circuit();
        qc.reset_circuit();
        qc << aqs::Y{i};
        qc.generate_circuit();
        af::array pauli_y = qc.circuit();
        qc.reset_circuit();
        qc << aqs::Z{i};
        qc.generate_circuit();
        af::array pauli_z = qc.circuit();
        qc.reset_circuit();

        coeff.push_back(0.25f * hilsch_product(identity, hamiltonian));
        coeff.push_back(0.25f * hilsch_product(pauli_x, hamiltonian));
        coeff.push_back(0.25f * hilsch_product(pauli_y, hamiltonian));
        coeff.push_back(0.25f * hilsch_product(pauli_z, hamiltonian));
    }

    return coeff;
}

std::vector<std::pair<std::string, af::cfloat>> decompose(uint32_t qubits, const af::array& hamiltonian)
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

        auto coeff = hs_product(temp, hamiltonian); 
        if (coeff != af::cfloat{0.f} || coeff != af::cfloat{-0.f})
            coeffs.push_back({str, coeff * (1.f / fast_pow2(qubits))});
    }

    return coeffs;
}

af::array compose(uint32_t qubits, const std::vector<std::pair<std::string, af::cfloat>>& data)
{
    af::array out = af::constant(af::cfloat{}, fast_pow2(qubits), fast_pow2(qubits));
    static af::cfloat ivals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {1.f,0.f} };
    static af::cfloat xvals[] = { {0.f,0.f}, {1.f,0.f}, {1.f,0.f}, {0.f,0.f} };
    static af::cfloat yvals[] = { {0.f,0.f}, {0.f,1.f}, {0.f,-1.f}, {0.f,0.f} };
    static af::cfloat zvals[] = { {1.f,0.f}, {0.f,0.f}, {0.f,0.f}, {-1.f,0.f} };

    static af::array I(2, 2, ivals);
    static af::array X(2, 2, xvals);
    static af::array Y(2, 2, yvals);
    static af::array Z(2, 2, zvals);
    for (const auto& pair : data)
    {
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

void add_matrix(aqs::QCircuit& qc, const std::string& shape, float coeff)
{
    if (qc.qubit_count() != shape.length())
        throw std::invalid_argument{"Qubit count and shape length must match"};

    if (shape.find_first_of("xyz") == std::string::npos)
        return;

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

    qc << aqs::RotZ{previous, coeff};

    //for (auto it = shape.rbegin() + (shape.length() - 1 - previous); it !)

    auto current = shape.find_last_of("xyz", previous - 1);
    while (current != std::string::npos && current != previous)
    {
        qc << aqs::CX(current, previous);
        previous = current;
        current = shape.find_last_of("xyz", previous);
    }
}

aqs::QCircuit hamiltonian_circuit(const af::array& hamiltonian, uint32_t steps)
{
    uint32_t qubits = fast_log2(hamiltonian.dims()[0]);
    auto hamiltonian_decomposition = decompose(qubits, hamiltonian);

    aqs::QCircuit qc{qubits};

    for (const auto& pair : hamiltonian_decomposition)
    {
        const std::string& shape = pair.first;
        af::cfloat coeff = pair.second;

        if (coeff.imag != 0.f)
            throw std::invalid_argument{"Must have a real hamiltonian"};

        add_matrix(qc, shape, coeff.real / (float) steps);
        qc << aqs::Barrier{false};
    }
    qc.generate_circuit();
    return qc;
}

void vqe_example()
{
    const uint32_t qubits = 2;
    const uint32_t depth = qubits;

    af::cfloat hamil_values[] = {
        {1.f},{2.f},{3.f},{4.f}
    };

    struct Data_t {
        aqs::QCircuit hamiltonian;
    };

    af::array hamil_matrix = af::diag(af::array(4, hamil_values), 0, false);

    auto hamil_circuit = hamiltonian_circuit(hamil_matrix, 1);
    Data_t data{hamil_circuit};
    //data.hamiltonian = hamil_circuit;

    std::vector<double> state_params(qubits * depth * 2, 1.);

    auto cost_function = [](const std::vector<double>& x, std::vector<double>& gradient, void* data) {
        static std::vector<float> buff(x.size());
        Data_t d = *static_cast<Data_t*>(data);
        const aqs::QCircuit& hamiltonian = d.hamiltonian;
        uint32_t qubits = hamiltonian.qubit_count();

        for (std::size_t i = 0; i < x.size(); ++i)
            buff[i] = static_cast<float>(x[i]);

        aqs::QCircuit qc{qubits};
        aqs::QSimulator qs(qubits);
        qc << aqs::Gate{linear_vqe(qubits, qubits, buff), 0};
        qc.generate_circuit();
        qs.generate_global_state();
        qs.simulate(qc);

        auto bra_state = af::transpose(qs.global_state(), true);

        qc << aqs::Gate{hamiltonian, 0};
        qc.generate_circuit();
        qs.generate_global_state();
        qs.simulate(qc);

        auto ket_state = qs.global_state();

        af::cfloat expectation = af::matmul(bra_state, ket_state)(0).scalar<af::cfloat>();

        double angle = std::acos(expectation.real);
        //if (expectation.imag < 0.f)
            //angle = -angle;

        return angle;
    };

    nlopt::opt opt(nlopt::LN_COBYLA, state_params.size());
    opt.set_min_objective(cost_function, &data);
    opt.set_xtol_rel(1e-3);
    opt.set_ftol_rel(1e-3);
    opt.set_maxeval(10000);
    opt.set_maxtime(20);

    double result = 0.0;
    std::cout << state_params;
    opt.optimize(state_params, result);

    std::vector<float> fparams(state_params.size());
    for (std::size_t i = 0; i < state_params.size(); ++i)
        fparams[i] = static_cast<float>(state_params[i]);

    auto state_circuit = linear_vqe(qubits, qubits, fparams);
    aqs::QSimulator qs(qubits);
    qs.simulate(state_circuit);

    std::cout << result << "\n";
    std::cout << state_params;
    aqs::print_global_state(qs);

    qs.simulate(hamil_circuit);
    aqs::print_global_state(qs);

    std::cout << "--------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    begin();
    optimize();

    vqe_example();

    uint32_t qubits = 4;
    uint32_t param_count = qubits * qubits * 2;
    auto h = hamiltonian(qubits);

    Data data;
    data.qubits = qubits;

    auto min_expectation = [](const std::vector<double>& x, std::vector<double>& gradient, void* data) {
        static std::vector<float> temp(x.size());
        const Data& d = *static_cast<Data*>(data);
        for (std::size_t i = 0; i < x.size(); ++i)
            temp[i] = static_cast<float>(x[i]);

        aqs::QCircuit qc{d.qubits};
        aqs::QSimulator qs(d.qubits);
        qc << aqs::Gate{linear_vqe(d.qubits, d.qubits, temp), 0};
        qc.generate_circuit();
        qs.generate_global_state();
        qs.simulate(qc);

        auto bra_state = af::transpose(qs.global_state(), true);

        qc << aqs::Gate{hamiltonian(d.qubits), 0};
        qc.generate_circuit();
        qs.generate_global_state();
        qs.simulate(qc);

        auto ket_state = qs.global_state();

        double expectation = static_cast<double>(af::real(af::matmul(bra_state, ket_state))(0).scalar<float>());

        return expectation;
    };

    std::vector<double> params(param_count);
    for (auto& val : params)
        val = dist(rd);

    nlopt::opt opt{nlopt::LN_COBYLA, param_count};
    opt.set_min_objective(min_expectation, &data);
    opt.set_xtol_rel(1e-3);
    opt.set_ftol_rel(1e-3);
    opt.set_maxeval(10000);
    opt.set_maxtime(20);
    //opt.set_initial_step(1e-2);
    //opt.set_lower_bounds(-af::Pi);
    //opt.set_upper_bounds(af::Pi);

    double result = 0.0;
    std::cout << params;
    //opt.optimize(params, result);

    std::vector<float> fparams(params.size());
    for (std::size_t i = 0; i < params.size(); ++i)
        fparams[i] = static_cast<float>(params[i]);

    auto state_circuit = linear_vqe(qubits, qubits, fparams);
    aqs::QSimulator qs(qubits);
    qs.simulate(state_circuit);

    std::cout << result << "\n";
    std::cout << params;
    aqs::print_global_state(qs);

    state_circuit << aqs::Gate(hamiltonian(qubits), 0);
    state_circuit.generate_circuit();
    qs.generate_global_state();
    qs.simulate(state_circuit);
    aqs::print_global_state(qs);

    af::cfloat mat_vals[] = {
        //{ 0.f } , {0.f} , {0.f} , {1.f}, { 0.f},{0.f,1.f},{0.f,0.f},{0.f}
        {1.f},{1.f},{1.f},{1.f}, {0.f,1.f},{0.f, 1.f},{0.f, 1.f},{0.f,1.f}
    };
    af::cfloat other[] = {
        {0.f}, {0.f}, {0.f}, {0.f, -1.f},
        {0.f}, {0.f}, {0.f,1.f}, {0.f},
        {0.f}, {0.f,-1.f}, {0.f}, {0.f},
        {0.f,1.f}, {0.f}, {0.f}, {0.f},
    };

    af::array mat = af::diag(af::array(8, mat_vals), 0, false);
    af_print(mat);
    //std::cout << get_coefficients(2, mat);
    auto d = decompose(3, mat);
    std::cout << d << std::endl;
    af_print(compose(3, d));

    af::array othermat = af::array(4, 4, other).T() * 2.f;
    af_print(mat);
    auto otherd = decompose(2, othermat);
    std::cout << otherd << std::endl;
    af_print(compose(2, otherd));


    af::cfloat ham_mat[] = {
        {0.f}, {0.f}, {0.f}, {1.f}
    };

    auto basic_hamiltonian = af::diag(af::array(4, ham_mat), 0, false);
    std::cout << decompose(2, basic_hamiltonian);
    auto hamil_circuit = hamiltonian_circuit(basic_hamiltonian, 1000);
    hamil_circuit.generate_circuit();
    aqs::print_circuit_text_image(hamil_circuit, aqs::QSimulator{2});
    aqs::print_circuit_matrix(hamil_circuit);
    
    aqs::QCircuit tempqc(2);
    for (int i = 0; i < 10000; ++i)
        tempqc << aqs::Gate(hamil_circuit, 0);

    tempqc.generate_circuit();
    aqs::print_circuit_matrix(tempqc);

    return 0;
}