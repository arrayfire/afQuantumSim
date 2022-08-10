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

#include "quantum.h"
#include "quantum_algo.h"
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
        stream << "{ " << val.first << " , " << val.second << " } , ";
    stream << "\b\b]\n";
    return stream;
}

void generic_hamiltonian_example()
{
    std::cout << "***** VQE example for a Hamiltonian *****\n\n";
    const uint32_t qubits = 2;
    const uint32_t depth = qubits;

    // Create a diagonal matrix as the hamiltonian
    af::cfloat hamil_values[] = {
        {0.f},{0.f},{0.f},{-1.f}
    };

    af::array hamil_matrix = af::diag(af::array(4, hamil_values), 0, false);

    std::cout << "Hamiltonian to find minimum eigenvalue of:\n";
    af_print(hamil_matrix);

    // Create a hamiltonian evolution of 1/10 the scale
    int scale = 10;
    auto hamil_circuit = aqs::hamiltonian_evolution_circuit(hamil_matrix, scale, true);
    std::cout << "Pauli representation of the hamiltonian: " << aqs::decompose_hamiltonian(hamil_matrix, 2) << "\n";

    std::cout << "Hamiltonian evolution circuit matrix:\n";
    aqs::print_circuit_matrix(hamil_circuit);
    
    std::cout << "\nHamiltonian evolution circuit representation:\n";    
    aqs::print_circuit_text_image(hamil_circuit, aqs::QSimulator{2});

    // Use variational quantum eigensolver to search for the smallest eigenvalue in the range [-10,10]
    auto pair = aqs::variational_quantum_eigensolver(hamil_matrix, (float)scale, aqs::VQE::LINEAR);
    const auto& result = pair.first;
    const auto& fparams = pair.second;

    // Display the state generator circuit from the paramaters
    auto state_circuit = aqs::linear_entanglement_varstate(qubits, qubits, fparams, true);
    aqs::QSimulator qs(qubits);
    qs.simulate(state_circuit);

    std::cout << "\nMinimum eigenvalue: " << result << "\n";
    std::cout << "\nAngle parameters" << fparams << "\n";
    std::cout << "\nEigenstate:\n";
    aqs::print_statevector(qs);

    qs.simulate(hamil_circuit);
    std::cout << "Resulting state:\n";
    aqs::print_statevector(qs);

    std::cout << "\n--------------------\n\n";
}

void hydrogen_molecule_example()
{
    std::cout << "***** Hydrogen molecule ground energy ***** \n\n";

    // Create the pauli-decomposition for the hydrogen molecule hamiltonian
    std::vector<std::pair<std::string, af::cfloat>> hamiltonian_description = {
        {"iizi", {-0.24274501250395486f}},
        {"iiiz", {-0.24274501250395486f}},
        {"iiii", {-0.04207255204090424f}},
        {"ziii", {0.17771358235540047f}},
        {"izii", {0.1777135823554005f}},
        {"zizi", {0.12293330446049033f}},
        {"iziz", {0.12293330446049033f}},
        {"ziiz", {0.16768338851167847f}},
        {"izzi", {0.16768338851167847f}},
        {"zzii", {0.17059759275420894f}},
        {"iizz", {0.1762766138632343}},
        {"yyxx", {-0.044750084051188126f}},
        {"xxyy", {-0.044750084051188126f}},
        {"yxxy", {0.044750084051188126f}},
        {"xyyx", {0.044750084051188126f}}
    };

    const uint32_t qubits = 4;
    const uint32_t depth = qubits;
    const int scale = 10;

    // Build hamiltonian matrix from the decomposition
    af::array hamiltonian_matrix = aqs::compose_hamiltonian(hamiltonian_description);
    std::cout << "Hydrogen molecule decomposed hamiltonian:\n" << aqs::decompose_hamiltonian(hamiltonian_matrix, qubits);

    // Find the smallest eigenvalue (= ground energy) for the hamiltonian
    auto pair = aqs::variational_quantum_eigensolver(hamiltonian_matrix, scale, aqs::VQE::LINEAR, 0.0f, 1000);
    auto result = pair.first;
    auto& params = pair.second;

    auto circuit = aqs::hamiltonian_evolution_circuit(hamiltonian_matrix, scale, true);

    auto state_circuit = aqs::full_entanglement_varstate(qubits, qubits, params, true);
    aqs::QSimulator qs(qubits);
    qs.simulate(state_circuit);

    std::cout << "\nMininum eigenvalue: " << result << "\n";
    std::cout << "\nParameters:\n";
    std::cout << params << "\n";
    std::cout << "Eigenstate:\n";
    aqs::print_statevector(qs);

    std::cout << "\nResulting state:\n";
    qs.simulate(circuit);
    aqs::print_statevector(qs);

    std::cout << "--------------------\n\n";
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    generic_hamiltonian_example();

    hydrogen_molecule_example();

    return 0;
}