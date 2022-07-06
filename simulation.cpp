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
    auto measures = qs.probabilities();

    double cost = 0.0;
    for (int i = 0; i < 3; ++i)
        cost += std::abs(measures[i] - target[i]);

    return cost;
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
    opt.set_stopval(1e-6);
    opt.set_xtol_abs(1e-5);
    opt.set_maxeval(10000);
    opt.set_maxtime(20);
    opt.set_initial_step(1e-1);

    double pi = af::Pi;
    opt.set_lower_bounds({-pi, -pi, -pi});
    opt.set_upper_bounds({pi, pi, pi});

    std::uniform_real_distribution<double> d{-pi, pi};
    std::vector<double> out = { d(rd) , d(rd), d(rd) };
    out = opt.optimize(out);

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
            qc << aqs::RotY{i, values[i]};
        for (uint32_t i = 0; i < qubits; ++i)
            qc << aqs::RotZ{i, values[qubits + i]};
    }

    qc.generate_circuit();
    return qc;
}

aqs::QCircuit hamiltonian(uint32_t qubits)
{
    aqs::QCircuit qc{qubits};

    for (uint32_t i = 0; i < qubits; ++i)
        qc << aqs::Z{i};

    qc.generate_circuit();
    return qc;
}

struct Data
{
    uint32_t qubits;
};

int main(int argc, char** argv)
{
    //begin();
    //optimize();

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
    opt.set_stopval(1e-6);
    opt.set_xtol_abs(1e-5);
    opt.set_maxeval(10000);
    opt.set_maxtime(20);
    opt.set_initial_step(1e-1);
    opt.set_lower_bounds(-af::Pi);
    opt.set_upper_bounds(af::Pi);

    double result = 0.0;
    opt.optimize(params, result);

    std::vector<float> fparams(params.size());
    for (std::size_t i = 0; i < params.size(); ++i)
        fparams[i] = static_cast<float>(params[i]);

    auto state_circuit = linear_vqe(qubits, qubits, fparams);
    aqs::QSimulator qs(qubits);
    qs.simulate(state_circuit);

    std::cout << result << "\n";
    std::cout << params;
    aqs::print_global_state(qs);

    return 0;
}