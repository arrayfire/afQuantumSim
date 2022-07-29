/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum_visuals.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_set>

namespace aqs
{

void print_state(const QState& state)
{
    std::cout << std::setprecision(4) << std::fixed << std::showpos;
    std::cout << state[0] << " |0> + " << state[1] << " |1>\n";
    std::cout << std::setprecision(7) << std::defaultfloat << std::noshowpos;
}

void print_statevector(const QSimulator& simulator)
{
    const int qubits = simulator.qubit_count();
    const int states = simulator.state_count();
    std::vector<af::cfloat> vals(states);

    simulator.statevector().host(vals.data());
    std::cout << std::setprecision(3) << std::fixed << std::showpos;
    for (int i = 0; i < states - 1; ++i)
        std::cout << vals[i] << "|" << binary_string(i, qubits) << "> + " << ((i + 1) % 4 ? "" : "\n");
    std::cout << vals[states - 1] << "|" << binary_string(states - 1, qubits) << ">\n\n";
    std::cout << std::setprecision(7) << std::defaultfloat << std::noshowpos;
}

void print_circuit_matrix(const QCircuit& circuit)
{
    std::cout << std::setprecision(3) << std::fixed << std::showpos;
    af::print("Circuit: ", circuit.circuit());
    std::cout << std::setprecision(7) << std::defaultfloat << std::noshowpos;
}

void print_profile(const std::array<uint32_t, 2>& profile)
{
    int rep_count = profile[0] + profile[1];
    std::cout << std::setprecision(3) << std::fixed << std::noshowpos;
    std::cout << "|0>: " << profile[0] * 100.f / rep_count << "% (" << profile[0] << ")\n" 
              << "|1>: " << profile[1] * 100.f / rep_count << "% (" << profile[1] << ")\n";
    std::cout << std::setprecision(7) << std::defaultfloat << std::noshowpos;
}

void print_profile(const std::vector<uint32_t>& profile)
{
    int qubits = fast_log2(profile.size());
    int rep_count = std::accumulate(profile.begin(), profile.end(), 0);
    std::cout << std::setprecision(2) << std::fixed;
    for (int i = 0; i < profile.size(); i++)
        std::cout << "|" << binary_string(i, qubits) << ">: " << std::setw(5) << profile[i] * 100.f / rep_count << "% (" << profile[i] << ")\n";
    std::cout << std::setprecision(7) << std::defaultfloat;
}

std::string parse_circuit_representation(std::string parse_string)
{
    std::regex info_regex{ "([0-9]+);(([0-9]+,[0-9];)+)" };
    std::regex initial_states_regex{ "([0-9]+),([0|1]);" };
    //std::regex gate_regex{ "([A-Za-z0-9]+),([0-9]+),([0-9]+):(([0-9]+,?)+);" }; // ASCII ONLY
    std::regex gate_regex{ "([^,|;]+),([0-9]+),([0-9]+):(([0-9]+,?)+);" }; // UTF8
    std::regex barrier_regex{ "[B|P];" };
    std::regex qubit_list_regex{ "([0-9]+),?" };

    // Remove whitespace
    parse_string.erase(std::remove_if(parse_string.begin(), parse_string.end(),
                       [](unsigned char c){
                           return c == '\n' || c == '\r' || c == ' ' || c == '\t';
                       }),
                       parse_string.end());

    uint32_t qubit_count = 0;
    std::array<bool, aqs::max_qubit_count> states{};

    // Buffer used to store the list of all control qubits in a gate
    std::vector<uint32_t> control_qubit_buffer;
    control_qubit_buffer.reserve(aqs::max_qubit_count);

    // Buffer used to store the list of all target qubits in a gate
    std::vector<uint32_t> target_qubit_buffer;
    target_qubit_buffer.reserve(aqs::max_qubit_count);

    // Set used to track the qubits in a list and assure uniqueness
    std::unordered_set<uint32_t> qubit_set;
    qubit_set.reserve(qubit_count);

    std::smatch info_match{};
    if (std::regex_search(parse_string, info_match, info_regex))
    {
        qubit_count = static_cast<uint32_t>(std::stoul(info_match[1].str()));

        if (!qubit_count)
            throw std::out_of_range{ "Circuit must contain at least one qubit" };

        if (qubit_count > aqs::max_qubit_count)
            throw std::out_of_range{ "Circuit exceeds supported max qubit count" };

        std::string initial_states_str = info_match[2].str();

        auto& initialized_qubits = qubit_set;

        std::smatch inital_states_match{};
        while (std::regex_search(initial_states_str, inital_states_match, initial_states_regex))
        {
            // Find the index of qubit state being set
            uint32_t qubit_index = static_cast<uint32_t>(std::stoul(inital_states_match[1].str()));

            if(initialized_qubits.count(qubit_index))
                throw std::invalid_argument{ "Cannot set the state of the same qubit multiple times" };

            initialized_qubits.insert(qubit_index);

            // Store the state the qubit is in
            auto state = static_cast<bool>(std::stoul(inital_states_match[2].str()));

            if (!(state == 1 || state == 0))
                throw std::invalid_argument{ "Invalid initial state" };

            states[qubit_index] = state;

            initial_states_str = inital_states_match.suffix();
        }

        initialized_qubits.clear();
    }

    parse_string = info_match.suffix();

    std::vector<std::stringstream> lines{};
    std::vector<std::size_t> lines_length{};

    const std::size_t line_count = qubit_count * 3;
    lines.resize(line_count);
    lines_length.resize(qubit_count, 0);

    for (uint32_t i = 0; i < qubit_count; ++i)
    {
        lines[i * 3    ] << "   ";
        lines[i * 3 + 1] << "|" << (int)states[i] << "⟩";
        lines[i * 3 + 2] << "   ";

        lines_length[i] += 3;
    }

    std::size_t current_circuit_length = 3;

    //Fills all the lines until all have the same width, adds tails if defined
    auto fill_circuit =
        [&lines, &lines_length, qubit_count, &current_circuit_length]
        (std::string center = "", std::string sides = "", std::size_t tail_len = 0)
    {
        for (uint32_t i = 0; i < qubit_count; ++i)
        {
            lines[i * 3    ] << repeat(current_circuit_length - lines_length[i], " ") << sides;
            lines[i * 3 + 1] << repeat(current_circuit_length - lines_length[i], "─") << center;
            lines[i * 3 + 2] << repeat(current_circuit_length - lines_length[i], " ") << sides;

            lines_length[i] = current_circuit_length + tail_len;
        }

        current_circuit_length += tail_len;
    };

    //Fills all the lines in the range until all have the same width
    auto fill_qubit_range =
        [&lines, &lines_length, qubit_count, &current_circuit_length]
        (uint32_t first_qubit, uint32_t last_qubit)
    {
        std::size_t max_line_length = 0;

        for (std::size_t i = first_qubit; i <= last_qubit; ++i)
            max_line_length = (lines_length[i] > max_line_length ? lines_length[i] : max_line_length);

        for (std::size_t i = first_qubit; i <= last_qubit; ++i)
        {
            auto line = i * 3;
            lines[line    ] << repeat(max_line_length - lines_length[i], " ");
            lines[line + 1] << repeat(max_line_length - lines_length[i], "─");
            lines[line + 2] << repeat(max_line_length - lines_length[i], " ");

            lines_length[i] = max_line_length;
        }

        current_circuit_length = (current_circuit_length < max_line_length ? max_line_length : current_circuit_length);
    };

    while (!parse_string.empty())
    {
        std::smatch gate_match{};

        if (std::regex_search(parse_string, gate_match, barrier_regex) && gate_match.prefix().str().empty())
        {
            if (gate_match.str() == "B;")
                fill_circuit("▒─", "▒ ", 2);
            else if (gate_match.str() == "P;")
                fill_circuit();
            else
                throw std::runtime_error{ "Unexpected error" };
        }
        else if (std::regex_search(parse_string, gate_match, gate_regex))
        {
            std::string gate_name = gate_match[1];
            uint32_t control_qubit_count = static_cast<uint32_t>(std::stoul(gate_match[2]));
            uint32_t target_qubit_count = static_cast<uint32_t>(std::stoul(gate_match[3]));
            std::string qubit_list_str = gate_match[4];
            std::smatch qubit_list_match{};

            while (std::regex_search(qubit_list_str, qubit_list_match, qubit_list_regex))
            {
                uint32_t qubit = static_cast<uint32_t>(std::stoul(qubit_list_match[1]));

                if (qubit >= qubit_count)
                    throw std::out_of_range{ "Qubit index is outside the circuit dimensions" };

                if (qubit_set.count(qubit))
                    throw std::invalid_argument{ "Invalid gate registers: cannot have the same qubit as multiple gate registers" };

                qubit_set.insert(qubit);

                if (control_qubit_buffer.size() < control_qubit_count)
                    control_qubit_buffer.push_back(qubit);
                else
                    target_qubit_buffer.push_back(qubit);

                qubit_list_str = qubit_list_match.suffix();
            }

            // Checking qubit lists
            if (control_qubit_buffer.size() != control_qubit_count ||
                target_qubit_buffer.size() != target_qubit_count)
                throw std::invalid_argument{ "Number of control and target qubits must match gate declaration" };

            // General Constraints
            if (!target_qubit_buffer.size())
                throw std::invalid_argument{ "Gate must contain at least one target qubit" };

            // Special Gate Contraints
            if (gate_name == "Swap" && target_qubit_buffer.size() != 2)
                throw std::invalid_argument{ "Number of target qubits for a Swap gate must be 2" };

            std::sort(control_qubit_buffer.begin(), control_qubit_buffer.end());
            std::sort(target_qubit_buffer.begin(), target_qubit_buffer.end());

            bool is_control_first = false;
            bool is_control_last = false;

            const uint32_t first_qubit = [&](){
                if (!control_qubit_buffer.empty())
                {
                    is_control_first = target_qubit_buffer.front() > control_qubit_buffer.front();
                    return is_control_first ? control_qubit_buffer.front() : target_qubit_buffer.front();
                }
                else
                {
                    is_control_first = false;
                    return target_qubit_buffer.front();
                }
            }();

            const uint32_t last_qubit = [&](){
                if (!control_qubit_buffer.empty())
                {
                    is_control_last = target_qubit_buffer.back() < control_qubit_buffer.back();
                    return is_control_last ? control_qubit_buffer.back() : target_qubit_buffer.back();
                }
                else
                {
                    is_control_last = false;
                    return target_qubit_buffer.back();
                }
            }();

            const std::string pad = " ";
            const std::string gate_str = pad + gate_name + pad;
            const std::size_t strlen = utf8str_len(gate_str);
            const std::size_t strmid = (strlen + 1) / 2 - 1;

            std::size_t len = 1;
            std::size_t mid = 1;

            if (gate_name == "Swap")
            {
                len = 3;
                mid = 1;
            }
            else
            {
                len = strlen + 6;
                mid = strmid + 3;
            }

            // Pad the lines in the range of the gate
            fill_qubit_range(first_qubit, last_qubit);

            uint32_t current_ctrl_index = control_qubit_buffer.empty() ? -1 : 0;
            uint32_t current_trgt_index = 0;

            uint32_t current_qubit = first_qubit;
            uint32_t previous_qubit = first_qubit;

            // Add marks until all target and control marks have been added
            while (current_ctrl_index < control_qubit_buffer.size() || current_trgt_index < target_qubit_buffer.size())
            {
                // Find the next current qubit in descending order
                if (current_ctrl_index < control_qubit_buffer.size())
                {
                    if (current_trgt_index < target_qubit_buffer.size())
                    {
                        current_qubit = target_qubit_buffer[current_trgt_index] < control_qubit_buffer[current_ctrl_index] ?
                                        target_qubit_buffer[current_trgt_index] : control_qubit_buffer[current_ctrl_index];
                    }
                    else
                    {
                        current_qubit = control_qubit_buffer[current_ctrl_index];
                    }
                }
                else
                {
                    current_qubit = target_qubit_buffer[current_trgt_index];   
                }

                // Fill sections in between the control qubits and gates
                for (uint32_t i = previous_qubit + 1; i < current_qubit; ++i)
                {
                    lines[i * 3    ] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");
                    lines[i * 3 + 1] << repeat(mid, "─") << "┼" << repeat(len - mid - 1, "─");
                    lines[i * 3 + 2] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");

                    lines_length[i] += len;
                }

                // Handle adding marks for the target qubit
                if (current_trgt_index < target_qubit_buffer.size() &&
                    target_qubit_buffer[current_trgt_index] == current_qubit)
                {
                    // Special case for SWAP gate
                    if (gate_name == "Swap")
                    {
                        uint32_t top_swap = target_qubit_buffer[0];
                        uint32_t bottom_swap = target_qubit_buffer[1];

                        if (top_swap > bottom_swap)
                            std::swap(top_swap, bottom_swap);

                        if (top_swap == current_qubit)
                        {
                            // Handle cases with control qubits
                            if (is_control_first)
                                lines[top_swap * 3] << " │ ";
                            else
                                lines[top_swap * 3] << "   ";

                            lines[top_swap * 3 + 1] << "─╳─";
                            lines[top_swap * 3 + 2] << " │ ";
                        }
                        else if (bottom_swap == current_qubit)
                        {
                            lines[bottom_swap * 3    ] << " │ ";
                            lines[bottom_swap * 3 + 1] << "─╳─";

                            // Handle cases with control qubits
                            if (is_control_last)
                                lines[bottom_swap * 3 + 2] << " │ ";
                            else
                                lines[bottom_swap * 3 + 2] << "   ";
                        }
                        else
                            throw std::runtime_error{ "Unexpected error" };

                        ++current_trgt_index;
                        lines_length[current_qubit] += len;
                        previous_qubit = current_qubit;
                    }
                    // Handle general gates
                    else
                    {
                        const uint32_t begin = current_qubit;
                        uint32_t count = 1;

                        // Find if there are multiple multi-qubit target qubit gates
                        for (uint32_t i = current_trgt_index + 1;
                            i < target_qubit_buffer.size() && target_qubit_buffer[i] == target_qubit_buffer[i - 1] + 1;
                            ++i)
                                ++count;

                        current_trgt_index += count;

                        // Top qubit marks
                        if (!is_control_first && current_qubit == target_qubit_buffer.front())
                            lines[begin * 3] << "  ┌" << repeat(strlen, "─") << "┐  ";
                        else
                            lines[begin * 3] << "  ┌" << repeat(strmid, "─") << "┴"
                                                << repeat(strlen - strmid - 1, "─") << "┐  ";

                        lines[begin * 3 + 1] << "──┤" << gate_str << "├──";

                        // Mid qubit marks
                        for (uint32_t i = begin + 1; i < begin + count - 1; ++i)
                        {
                            lines[i * 3    ] << "  │" << repeat(strlen, " ") << "│  ";
                            lines[i * 3 + 1] << "──┤" << repeat(strlen, " ") << "├──";
                            lines[i * 3 + 2] << "  │" << repeat(strlen, " ") << "│  ";
                        }

                        // Bottom qubit marks
                        if (count > 1)
                        {
                            lines[(begin            ) * 3 + 2] << "  │" << repeat(strlen, " ") << "│  ";
                            lines[(begin + count - 1) * 3    ] << "  │" << repeat(strlen, " ") << "│  ";
                            lines[(begin + count - 1) * 3 + 1] << "──┤" << repeat(strlen, " ") << "├──";
                        }
                        if (!is_control_last && current_qubit + count - 1 == target_qubit_buffer.back())
                            lines[(begin + count - 1) * 3 + 2] << "  └" << repeat(strlen, "─") << "┘  ";
                        else
                            lines[(begin + count - 1) * 3 + 2] << "  └" << repeat(strmid, "─") << "┬"
                                                                << repeat(strlen - strmid - 1, "─") << "┘  ";

                        for (uint32_t i = begin; i < begin + count; ++i)
                            lines_length[i] += len;

                        previous_qubit = begin + count - 1;
                    }
                }
                else if (current_ctrl_index < control_qubit_buffer.size() &&
                         control_qubit_buffer[current_ctrl_index] == current_qubit)
                {
                    // If control qubit is the first qubit in the gate
                    if (is_control_first && current_qubit == control_qubit_buffer.front())
                    {
                        lines[current_qubit * 3.   ] << repeat(len, " ");
                        lines[current_qubit * 3 + 1] << repeat(mid, "─") << "█" << repeat(len - mid - 1, "─");
                        lines[current_qubit * 3 + 2] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");
                    }
                    // If control qubit is the last qubit in the gate
                    else if (is_control_last && current_qubit == control_qubit_buffer.back())
                    {
                        lines[current_qubit * 3    ] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");
                        lines[current_qubit * 3 + 1] << repeat(mid, "─") << "█" << repeat(len - mid - 1, "─");
                        lines[current_qubit * 3 + 2] << repeat(len, " ");
                    }
                    // Control qubits in between
                    else
                    {
                        lines[current_qubit * 3    ] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");
                        lines[current_qubit * 3 + 1] << repeat(mid, "─") << "█" << repeat(len - mid - 1, "─");
                        lines[current_qubit * 3 + 2] << repeat(mid, " ") << "│" << repeat(len - mid - 1, " ");
                    }

                    lines_length[current_qubit] += len;
                    ++current_ctrl_index;

                    previous_qubit = current_qubit;
                }
                else
                    throw std::runtime_error{ "Unexpected error" };
            }

            // Update the circuit length if the line length is larger
            if (lines_length[first_qubit] > current_circuit_length)
                current_circuit_length = lines_length[first_qubit];

            // Reset the buffers for next iteration
            control_qubit_buffer.clear();
            target_qubit_buffer.clear();
            qubit_set.clear();
        }
        else
        {
            throw std::invalid_argument{ "Could not parse the given circuit string representation" };
        }

        // Continue parsing remaining gates
        parse_string = gate_match.suffix();
    }
    
    // Fill circuit with padding till end
    fill_circuit();

    //Add newlines
    std::size_t total_length = 0;
    for (auto& line : lines)
    {
        line << "\n";

        line.seekp(0, std::ios::end);
        total_length += line.tellp();
    }

    std::string output;
    output.reserve(total_length + 1);
    output.append(1, '\n');

    //Append all lines into output string
    for (const auto& line : lines)
        output.append(line.str());

    return output;
}

std::string gen_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator)
{
    if (circuit.qubit_count() != simulator.qubit_count())
        throw std::invalid_argument{ "Circuit and simulator qubit count must match" };

    //Add number of qubits
    const int qubits = circuit.qubit_count();

    std::stringstream circuit_representation;
    circuit_representation << qubits << ";";

    //Add states
    for (int i = 0; i < qubits; ++i)
    {
        const auto& q = simulator.qubit(i);
        int val = 0;
        if (q == aqs::QState::zero())
            val = 0;
        else if (q == aqs::QState::one())
            val = 1;
        else
            throw std::invalid_argument{"Superposed inital states not supported"};

        circuit_representation << i << "," << val << ";";
    }

    //Add circuit gates
    circuit_representation << circuit.representation();

    return parse_circuit_representation(circuit_representation.str());
}

std::string gen_circuit_text_image(std::string schematic)
{
    return parse_circuit_representation(std::move(schematic));
}

void print_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator)
{
    std::cout << gen_circuit_text_image(circuit, simulator) << std::endl;
}
}