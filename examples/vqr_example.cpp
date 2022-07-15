/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#include <iostream>
#include <vector>
#include <numeric>
#include <random>

#include <nlopt.hpp>

#include "quantum.h"
#include "quantum_algo.h"
#include "quantum_visuals.h"

std::string polynomial_string(uint32_t degree, const std::vector<double>& params)
{
    if (degree + 1 != params.size())
        throw std::invalid_argument{"Number of params must match the degree of the polynomial"};

    std::stringstream str;

    str << "y = ";
    if (params[0] != 0 && degree != 0)
        str << params[0];
    if (degree > 0 && params[1] != 0)
    {
        if (params[0] != 0)
            str << (params[1] > 0.0 ? " - " : " + ");
        if (std::abs(params[1]) == 1)
            str << "x";
        else
            str << std::abs(params[1]) << " x";
    }
    for (uint32_t i = 2; i < degree + 1; ++i)
    {
        if (params[i] == 1)
            str << " + " << "x ^ " << i;
        else if (params[i] == 0)
            continue;
        else
            str << (params[i] < 0.0 ? " - " : " + ") << std::abs(params[i]) << " x ^ " << i;
    }

    return str.str();
}

/**
 * @brief Finds the inner product of two real normalized vectors using quantum circuits
 * 
 * @note Length of the vectors must be equal and a power of two
 * 
 * @details Uses probability from qubit measurement to determine the inner product
 * 
 * @see [implementation source](https://qiskit.org/textbook/ch-demos/variational-quantum-regression.html)
 * 
 * @param lhs 
 * @param rhs 
 * @return double 
 */
double quantum_inner_product(const af::array& lhs, const af::array& rhs)
{
    if (lhs.dims()[0] != rhs.dims()[0] &&
        lhs.dims()[1] != 1 && rhs.dims()[1] != 1 &&
        lhs.dims()[2] != 1 && rhs.dims()[2] != 1 &&
        lhs.dims()[3] != 1 && rhs.dims()[3] != 1)
        throw std::invalid_argument{"Cannot compute inner product between different size vectors"};

    if (af::norm(rhs) == 0.0)
        return 0.0;
        
    uint32_t size = lhs.dims()[0];
    uint32_t qubits = std::ceil(fast_log2(size));

    aqs::QCircuit qc{qubits + 1};
    qc << aqs::H{0};
    qc.generate_circuit();

    auto global_state = af::complex(af::join(0, lhs / af::norm(lhs), rhs / af::norm(rhs)));
    global_state /= af::norm(global_state);

    aqs::QSimulator qs{qubits + 1, global_state};
    qs.simulate(qc);

    auto result = qs.qubit_probability_false(0);

    return result * 2 - 1.0;
}

double linear_cost_function(const std::vector<double>& params, std::vector<double>& grad, void* data)
{
    auto a = params[1];
    auto b = params[0];

    const auto& input = static_cast<const af::array*>(data)[0];
    const auto& output = static_cast<const af::array*>(data)[1];

    auto output_norm = af::norm(output);

    // Predict output
    af::array predict = (input * a) + b;
    auto predict_norm = af::norm(predict);

    if (predict_norm == 0.0)
        return std::numeric_limits<double>::infinity();
    predict /= predict_norm;

    auto dir_diff = quantum_inner_product(output, predict);

    //return std::exp((1.0 - dir_diff) * (1.0 - dir_diff)) * (predict_norm / output_norm + output_norm / predict_norm);
    return (1.0 + (1.0 - dir_diff) * (1.0 - dir_diff)) * (predict_norm / output_norm + output_norm / predict_norm) / 2.0;
}

struct PolynomialData
{
    af::array& input;
    af::array& output;
    uint32_t degree;
};

/**
 * @brief Cost function for fitting a general nth degree polynomial
 * 
 * @details The minimum value for this function is 1.0
 *          The pointer pass to data should be a pointer to a PolynomialData variable
 * 
 * @param params 
 * @param grad 
 * @param data 
 * @return double 
 */
double polynomial_cost_function(const std::vector<double>& params, std::vector<double>& grad, void* data)
{
    const PolynomialData& temp = *static_cast<const PolynomialData*>(data);
    const auto& input = temp.input;
    const auto& output = temp.output;
    const auto& degree = temp.degree;

    auto output_norm = af::norm(output);

    //Determine the predicted output
    af::array predict = af::constant(params[0], input.dims()[0]);
    af::array current = input;
    for (uint32_t i = 1; i < degree + 1; ++i)
    {
        predict += current * params[i];
        current *= input;
    }

    auto predict_norm = af::norm(predict);
    if (predict_norm == 0.0)
        return std::numeric_limits<double>::infinity();

    predict /= predict_norm;

    auto dir_diff = quantum_inner_product(output, predict);

    //return std::exp((1.0 - dir_diff) * (1.0 - dir_diff)) * (predict_norm / output_norm + output_norm / predict_norm);
    return (1.0 + (1.0 - dir_diff) * (1.0 - dir_diff)) * (predict_norm / output_norm + output_norm / predict_norm) / 2.0;
}

void equality_regression()
{   
    std::cout << "**** Equality (y = x) example ****\n\n";

    std::random_device rd;
    af::randomEngine eng{AF_RANDOM_ENGINE_MERSENNE, rd()};
    int32_t points = 32;

    //Generate random points in the y=x line
    af::array vals[2];
    vals[0] = af::randu(points, f32, eng);
    vals[1] = vals[0];

    // Set optimization constraints
    nlopt::opt opt{nlopt::LN_NELDERMEAD, 2};
    opt.set_min_objective(linear_cost_function, &vals);
    opt.set_maxeval(200);

    double result = 0.0;
    std::vector<double> params(2, 0.5);

    // Initialize params randomly in range [-5.0 , 5.0]
    double range = 1.0;
    (af::randu(2, f64) * range * 2 - range).host(params.data());
    
    std::cout << "Expected line equation: " << polynomial_string(1, { 0.0 , 1.0 }) << "\n" <<
                 "Regression result " << polynomial_string(1, params) << "\n" <<
                 "Minimized cost function value: " << result << "\n\n";
    std::cout << "**************************\n\n";
}

void linear_regression()
{
    std::cout << "**** Linear example ****\n\n";

    std::random_device rd;
    af::randomEngine eng{AF_RANDOM_ENGINE_MERSENNE, rd()};

    //Generate random x points
    const int32_t points = 128;
    af::array vals[2];
    vals[0] = af::randu(points, f32, eng);

    //Evaluate the points in the line y
    const double a = 4.0;
    const double b = 5.0;
    vals[1] = (vals[0] * a) + b;

    // Set optimization constraints
    nlopt::opt opt{nlopt::LN_NELDERMEAD, 2};
    opt.set_min_objective(linear_cost_function, vals);
    opt.set_maxeval(1000);

    double result = 0.0;
    std::vector<double> params(2);

    // Initialize params randomly in range [-5.0 , 5.0]
    double range = 5.0;
    (af::randu(2, f64) * range * 2 - range).host(params.data());

    // Optimize and get param values and minimized cost function
    opt.optimize(params, result);

    std::cout << "Expected line equation: " << polynomial_string(1, { b , a }) << "\n" <<
                 "Regression result " << polynomial_string(1, params) << "\n" <<
                 "Minimized cost function value: " << result << "\n\n";
    std::cout << "**************************\n\n";
}

void polynomial_regression()
{
    std::cout << "**** Polynomial example ****\n\n";

    std::random_device rd;
    af::randomEngine eng{AF_RANDOM_ENGINE_MERSENNE, rd()};

    const uint32_t degree = 2;
    const double a = 4.0;
    const double b = -5.0;
    const double c = 3.5;

    //Generate random x points
    const uint32_t points = 64;
    af::array input = af::randu(points, f32, eng);

    //Generate the polynomial y values
    af::array output = (input * input) * a + input * b + c;

    PolynomialData data { input , output , degree };

    // Set optimization constraints
    nlopt::opt opt{nlopt::LN_NELDERMEAD, degree + 1};
    opt.set_min_objective(polynomial_cost_function, &data);
    opt.set_maxeval(1000);

    std::vector<double> params(degree + 1);
    double result = 0.0;

    // Initialize params randomly in range [-5.0 , 5.0]
    double range = 5.0;
    (af::randu(degree + 1, f64) * range * 2 - range).host(params.data());

    // Optimize and get param values and minimized cost function
    opt.optimize(params, result);

    std::cout << "Expected polynomial: " << polynomial_string(degree, { c , b , a }) << "\n" <<
                 "Regression result " << polynomial_string(degree, params) << "\n" <<
                 "Minimized cost function value: " << result << "\n\n";
    std::cout << "**************************\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << '\n';

    equality_regression();

    linear_regression();

    polynomial_regression();

    return 0;
}