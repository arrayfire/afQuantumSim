/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#pragma once

#include "quantum.h"

namespace aqs
{
    /**
     * @brief Prints the given qubit state to the stdout stream
     * 
     * @param state Qubit State
     */
    void print_state(const QState& state);

    /**
     * @brief Prints the circuit matrix represenatation to the stdout stream
     * 
     * @param circuit circuit to be displayed
     */
    void print_circuit_matrix(const QCircuit& circuit);

    /**
     * @brief Prints the global state stored in the simulator to the stdout stream
     * 
     * @param simulator simulator storing the global state
     */
    void print_statevector(const QSimulator& simulator);

    /**
     * @brief Prints the given qubit state profile to the stdout stream
     * 
     * @param profile profile of the measurements done on a qubit
     */
    void print_profile(const std::array<uint32_t, 2>& profile);

    /**
     * @brief Prints the given global state profile to the stdout stream
     * 
     * @param profile profile of the measurements done on the quantum simulation global state
     */
    void print_profile(const std::vector<uint32_t>& profile);

    /**
     * @brief Generates a text image using the circuit passed and initial states stored in simulator
     * 
     * @param circuit circuit to display
     * @param simulator simulator containing initial states
     * @return std::string utf8 text image
     */
    std::string gen_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator);

    /**
     * @brief Generates a text image using the schematic with the custom AQS displayer language
     * 
     * @param schematic string of the schematic
     * @return std::string utf8 text image
     */
    std::string gen_circuit_text_image(std::string schematic);

    /**
     * @brief Prints the text image of the circuit and inititial states passed to the stdout stream
     * 
     * @param circuit circuit to display
     * @param simulator simulator containing initial states
     */
    void print_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator);
}
