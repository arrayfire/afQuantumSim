/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h"
#include "quantum_gates.h"
#include "quantum_visuals.h"

#include <iostream>

void draw_schematic()
{
    std::string schematic = R"(
        3;

        0,1;
        1,0;

        H,0,1: 1;
        X,1,1: 1 , 0;
        Z,2,1: 0 , 2 , 1;
        Swap,0,2: 0 , 2;
    )";

    std::cout << aqs::gen_circuit_text_image(schematic) << std::endl;
}

void draw_circuit()
{
    uint32_t qubits = 3;

    // Create circuit
    aqs::QCircuit qc{ qubits };

    // Add gates
    qc << aqs::H{1}
       << aqs::CX{ 1 , 0 }
       << aqs::Gate(aqs::NControl_Gate(3, {0 , 2}, 1, aqs::Z::gate()), 0)
       << aqs::Swap{ 0 , 2 };

    // Create simulator
    aqs::QSimulator qs{ qubits , { aqs::QState::one() , aqs::QState::zero() , aqs::QState::zero() } };

    // Print circuit
    aqs::print_circuit_text_image(qc, qs);
}

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);

    draw_circuit();
    draw_schematic();

    return 0;
}
