/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/
#include "quantum_visuals.h"

#include <iostream>
#include <iomanip>
#include <numeric>

namespace aqs
{

void print_state(const QState& state)
{
    std::cout << std::setprecision(4) << std::fixed << std::showpos;
    std::cout << state[0] << " |0> + " << state[1] << " |1>\n";
    std::cout << std::setprecision(7) << std::defaultfloat << std::noshowpos;
}

void print_global_state(const QSimulator& simulator)
{
    const int qubits = simulator.qubit_count();
    const int states = simulator.state_count();
    std::vector<af::cfloat> vals(states);

    simulator.global_state().host(vals.data());
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

std::string parse_circuit_representation(const std::string& str)
{
    //First part: read the number of qubits of the circuit
    const auto count_end_pos =  str.find('\n', 0);
    if (count_end_pos == str.npos)
        throw std::invalid_argument{"string did not specify qubit count"};

    const auto qubit_count_str = str.substr(0, count_end_pos);

    const auto qubit_count = std::stoi(qubit_count_str);
    if (qubit_count < 1 || qubit_count > 31)
        throw std::invalid_argument{"invalid qubit quantity for circuit"};

    std::vector<std::string> lines;
    std::vector<std::size_t> lines_length;
    lines.resize(3 * qubit_count);
    lines_length.resize(3 * qubit_count, 0);

    //Second Part: Find the initial states that the circuit's qubits are intialized
    auto states_current_begin = count_end_pos + 1;
    auto states_current_end = states_current_begin;
    for (int i = 0; i < qubit_count; ++i)
    {
        const auto qubit_end = str.find(',', states_current_begin);
        states_current_end = str.find(';', states_current_begin);

        const auto qubit_str = str.substr(states_current_begin, qubit_end - states_current_begin);

        const auto qubit = std::stoi(qubit_str);
        if (qubit != i)
            throw std::invalid_argument{"invalid qubit ordering for qubit initialization"};

        const auto state_str = str.substr(qubit_end + 1, states_current_end - qubit_end - 1);

        const auto state = std::stoi(state_str);
        if (!(state == 0 || state == 1))
            throw std::invalid_argument{"invalid initial qubit state"};

        lines[i * 3    ].append("   ");
        lines[i * 3 + 1].append("|").append(std::to_string(state)).append("⟩");
        lines[i * 3 + 2].append("   ");
        lines_length[i * 3    ] += 3;
        lines_length[i * 3 + 1] += 3;
        lines_length[i * 3 + 2] += 3;

        states_current_begin = states_current_end + 1;
    }

    if (str.at(states_current_begin) != '\n')
        throw std::invalid_argument{"expected to contain newline separator at given position"};

    auto current_circuit_length = 3;
    auto gates_current_begin = states_current_begin + 1;
    auto gates_current_end = gates_current_begin;

    //Fills all the lines until all have the same width, adds tails if defined
    auto fill_circuit = [&lines, &lines_length, qubit_count, &current_circuit_length](std::string center = "", std::string sides = "", int tail_len = 0) {
        for (int i = 0; i < qubit_count; ++i)
        {
            lines[i * 3    ].append(current_circuit_length - lines_length[i * 3], ' ').append(sides);
            lines[i * 3 + 1].append(repeat(current_circuit_length - lines_length[i * 3 + 1], "─")).append(center);
            lines[i * 3 + 2].append(current_circuit_length - lines_length[i * 3 + 2], ' ').append(sides);
            lines_length[i * 3    ] = current_circuit_length + tail_len;
            lines_length[i * 3 + 1] = current_circuit_length + tail_len;
            lines_length[i * 3 + 2] = current_circuit_length + tail_len;
        }
        current_circuit_length += tail_len;
    };

    //Fills all the lines in the range until all have the same width
    auto fill_range = [&lines, &lines_length, qubit_count, &current_circuit_length](int begin, int end) {
        int max = 0;
        for (int i = begin; i < end; ++i)
            max = (lines_length[i * 3] > max ? lines_length[i * 3] : max);

        for (int i = begin; i < end; ++i)
        {
            lines[i * 3    ].append(max - lines_length[i * 3], ' ');
            lines[i * 3 + 1].append(repeat(max - lines_length[i * 3 + 1], "─"));
            lines[i * 3 + 2].append(max - lines_length[i * 3 + 2], ' ');
            lines_length[i * 3    ] = max;
            lines_length[i * 3 + 1] = max;
            lines_length[i * 3 + 2] = max;
        }
    };

    //Third Part: Parse and generate the circuit from the gates
    while(gates_current_begin != str.length())
    {
        auto gate_name_end = str.find(',', gates_current_begin);
        gates_current_end = str.find(';', gates_current_begin);
        if (gate_name_end > gates_current_end)
            gate_name_end = gates_current_end;
        std::string gate_name = str.substr(gates_current_begin, gate_name_end - gates_current_begin);
        auto control_qubit_begin = 0;
        auto control_qubit_count = 0;
        auto target_qubit_begin = 0;
        auto target_qubit_count = 0;
        
        //If a barrier is detected
        if (gate_name == "B")
        {
            fill_circuit("▒─", "▒ ", 2);
            gates_current_begin = gates_current_end + 1;
            continue;
        }
        //If padding is detected
        else if (gate_name == "P")
        {
            fill_circuit();
            gates_current_begin = gates_current_end + 1;
            continue;
        }
        //If any fundamental 1-qubit gate is detected
        else if (gate_name == "H" || gate_name == "X" || gate_name == "Y" || gate_name == "Z" || gate_name == "Phase" ||
                 gate_name == "T" || gate_name == "T†" || gate_name == "S" || gate_name == "S†")
        {
            const auto qubit_end = gates_current_end;
            const auto qubit_str = str.substr(gate_name_end + 1, qubit_end - gate_name_end);
            const auto qubit = std::stoi(qubit_str);
            target_qubit_begin = qubit;
            target_qubit_count = 1;
        }
        //If a Swap gate is detected
        else if (gate_name == "Swap")
        {
            const auto trgt_qubit1_end = str.find(',', gate_name_end + 1);
            const auto trgt_qubit1_str = str.substr(gate_name_end + 1, trgt_qubit1_end - gate_name_end);
            const auto trgt_qubit2_end = gates_current_end;
            const auto trgt_qubit2_str = str.substr(trgt_qubit1_end + 1, trgt_qubit2_end - trgt_qubit1_end);
            int qubit1 = std::stoi(trgt_qubit1_str);
            int qubit2 = std::stoi(trgt_qubit2_str);
            int top = qubit1 < qubit2 ? qubit1 : qubit2;
            int bottom = qubit1 < qubit2 ? qubit2 : qubit1;

            if (qubit1 > 31 || qubit2 > 31 || qubit1 < 0 || qubit2 < 0)
                throw std::invalid_argument{"Invalid target swap qubits"};

            fill_range(top, bottom + 1);

            lines[top * 3    ].append("   ");
            lines[top * 3 + 1].append("─╳─");
            lines[top * 3 + 2].append(" │ ");
            lines_length[top * 3    ] += 3;
            lines_length[top * 3 + 1] += 3;
            lines_length[top * 3 + 2] += 3;

            lines[bottom * 3    ].append(" │ ");
            lines[bottom * 3 + 1].append("─╳─");
            lines[bottom * 3 + 2].append("   ");
            lines_length[bottom * 3    ] += 3;
            lines_length[bottom * 3 + 1] += 3;
            lines_length[bottom * 3 + 2] += 3;
            for (int i = top + 1; i < bottom; ++i)
            {
                lines[i * 3    ].append(" │ ");
                lines[i * 3 + 1].append("─┼─");
                lines[i * 3 + 2].append(" │ ");
                lines_length[i * 3    ] += 3;
                lines_length[i * 3 + 1] += 3;
                lines_length[i * 3 + 2] += 3;
            }

            if (lines_length[top * 3] > current_circuit_length)
                current_circuit_length = lines_length[top * 3];

            gates_current_begin = gates_current_end + 1;
            continue;
        }
        //If a Control-Swap is detected
        else if (gate_name == "CSwap")
        {
            const auto ctrl_qubit_end = str.find(',', gate_name_end + 1);
            const auto ctrl_qubit_str = str.substr(gate_name_end + 1, ctrl_qubit_end - gate_name_end);
            const auto trgt_qubit1_end = str.find(',', ctrl_qubit_end + 1);
            const auto trgt_qubit1_str = str.substr(ctrl_qubit_end + 1, trgt_qubit1_end - ctrl_qubit_end);
            const auto trgt_qubit2_end = gates_current_end;
            const auto trgt_qubit2_str = str.substr(trgt_qubit1_end + 1, trgt_qubit2_end - trgt_qubit1_end);
            int ctrl = std::stoi(ctrl_qubit_str);
            int qubit1 = std::stoi(trgt_qubit1_str);
            int qubit2 = std::stoi(trgt_qubit2_str);
            int top = qubit1 < qubit2 ? qubit1 : qubit2;
            int bottom = qubit1 < qubit2 ? qubit2 : qubit1;

            if (ctrl < 0 || ctrl > 31)
                throw std::invalid_argument{"Invalid control qubit"};
            if (qubit1 > 31 || qubit2 > 31 || qubit1 < 0 || qubit2 < 0)
                throw std::invalid_argument{"Invalid target swap qubits"};

            if (ctrl < top)
            {
                fill_range(ctrl, bottom + 1);
                lines[top * 3    ].append(" │ ");
                lines[bottom * 3 + 2].append("   ");
                lines[ctrl * 3].append("   ");
                lines[ctrl * 3 + 2].append(" │ ");
                for (int i = ctrl + 1; i < top; ++i)
                {
                    lines[i * 3    ].append(" │ ");
                    lines[i * 3 + 1].append("─┼─");
                    lines[i * 3 + 2].append(" │ ");
                    lines_length[i * 3    ] += 3;
                    lines_length[i * 3 + 1] += 3;
                    lines_length[i * 3 + 2] += 3;
                }
            }
            else if (ctrl > bottom)
            {
                fill_range(top, ctrl + 1);
                lines[top * 3    ].append("   ");
                lines[bottom * 3 + 2].append(" │ ");
                lines[ctrl * 3].append(" │ ");
                lines[ctrl * 3 + 2].append("   ");
                for (int i = bottom + 1; i < ctrl; ++i)
                {
                    lines[i * 3    ].append(" │ ");
                    lines[i * 3 + 1].append("─┼─");
                    lines[i * 3 + 2].append(" │ ");
                    lines_length[i * 3    ] += 3;
                    lines_length[i * 3 + 1] += 3;
                    lines_length[i * 3 + 2] += 3;
                }
            }
            else
                throw std::invalid_argument{"Given Control qubit position for swap gate not supported"};

            lines[ctrl * 3 + 1].append("─█─");
            lines[top * 3 + 1].append("─╳─");
            lines[top * 3 + 2].append(" │ ");
            lines_length[top * 3    ] += 3;
            lines_length[top * 3 + 1] += 3;
            lines_length[top * 3 + 2] += 3;

            lines_length[ctrl * 3    ] += 3;
            lines_length[ctrl * 3 + 1] += 3;
            lines_length[ctrl * 3 + 2] += 3;

            lines[bottom * 3    ].append(" │ ");
            lines[bottom * 3 + 1].append("─╳─");
            lines_length[bottom * 3    ] += 3;
            lines_length[bottom * 3 + 1] += 3;
            lines_length[bottom * 3 + 2] += 3;
            for (int i = top + 1; i < bottom; ++i)
            {
                lines[i * 3    ].append(" │ ");
                lines[i * 3 + 1].append("─┼─");
                lines[i * 3 + 2].append(" │ ");
                lines_length[i * 3    ] += 3;
                lines_length[i * 3 + 1] += 3;
                lines_length[i * 3 + 2] += 3;
            }

            if (lines_length[top * 3] > current_circuit_length)
                current_circuit_length = lines_length[top * 3];

            gates_current_begin = gates_current_end + 1;
            continue;
        }
        //Any other gate is detected
        else
        {
            //Read the number of control qubits of the gate
            const auto control_qubit_count_end = str.find(',', gate_name_end + 1);
            const auto control_qubit_count_str = str.substr(gate_name_end + 1, control_qubit_count_end - gate_name_end + 1);
            control_qubit_count = std::stoi(control_qubit_count_str);

            //Read the number of target qubits of the gate
            const auto target_qubit_count_end = str.find(':', control_qubit_count_end + 1);
            const auto target_qubit_count_str = str.substr(control_qubit_count_end + 1, target_qubit_count_end - control_qubit_count_end - 1);
            target_qubit_count = std::stoi(target_qubit_count_str);
            
            if (target_qubit_count < 1)
                throw std::invalid_argument{"Custom gate must contain at least one target qubit"};
            if (target_qubit_count > 31)
                throw std::invalid_argument{"Invalid target qubit count"};
            if (control_qubit_count < 0 || target_qubit_count > 31)
                throw std::invalid_argument{"Invalid control qubit count"};

            auto ctrl_qubit_temp_begin = target_qubit_count_end + 1;
            auto ctrl_qubit_temp_end = str.find(',', ctrl_qubit_temp_begin);
            if (control_qubit_count > 0)
            {
                auto ctrl_qubit_temp_str = str.substr(ctrl_qubit_temp_begin, ctrl_qubit_temp_end - ctrl_qubit_temp_begin);
                control_qubit_begin = std::stoi(ctrl_qubit_temp_str);
                for (int i = 1; i < control_qubit_count; ++i)
                {
                    ctrl_qubit_temp_begin = ctrl_qubit_temp_end + 1;
                    ctrl_qubit_temp_end = str.find(',', ctrl_qubit_temp_begin);
                    auto ctrl_qubit = std::stoi(str.substr(ctrl_qubit_temp_begin, ctrl_qubit_temp_end - ctrl_qubit_temp_begin));
                    if (i + control_qubit_begin != ctrl_qubit)
                        throw std::invalid_argument{"Non-concurrent control qubits is not supported"};
                }
            }

            auto trgt_qubit_temp_begin = control_qubit_count != 0 ? ctrl_qubit_temp_end + 1 : ctrl_qubit_temp_begin;
            auto trgt_qubit_temp_end = control_qubit_count != 0 ? str.find(',', trgt_qubit_temp_begin) : ctrl_qubit_temp_end;
            target_qubit_begin = std::stoi(str.substr(trgt_qubit_temp_begin, trgt_qubit_temp_end - ctrl_qubit_temp_begin));
            for (int i = 1; i < target_qubit_count; ++i)
            {
                trgt_qubit_temp_begin = trgt_qubit_temp_end + 1;
                trgt_qubit_temp_end = str.find((i != target_qubit_count - 1) ? ',' : ';' , trgt_qubit_temp_begin);
                auto trgt_qubit = std::stoi(str.substr(trgt_qubit_temp_begin, trgt_qubit_temp_end - trgt_qubit_temp_begin));
                if (i + target_qubit_begin != trgt_qubit)
                    throw std::invalid_argument{"Non-concurrent target qubits is not supported"};
            }
        }

        const std::string padding = " ";

        //Remove C (control gate mark)
        for (int i = 0; i < control_qubit_count; ++i)
            if (gate_name[0] == 'C') gate_name.erase(0, 1);
        const std::string gate_str = padding + gate_name + padding;
        const std::size_t gate_strlen = utf8str_len(gate_str);
        const std::size_t gate_strmid = (gate_strlen + 1) / 2 - 1;
        const std::size_t offset = 3 + gate_strlen + 3;

        //Deal with circuits with control qubits normally
        if (control_qubit_count == 0)
        {
            fill_range(target_qubit_begin, target_qubit_begin + target_qubit_count);

            lines[target_qubit_begin * 3    ].append("  ┌").append(repeat(gate_strlen, "─")).append("┐  ");
            lines[target_qubit_begin * 3 + 1].append("──┤").append(gate_str).append("├──");

            lines_length[target_qubit_begin * 3    ] += offset;
            lines_length[target_qubit_begin * 3 + 1] += offset;
            if (target_qubit_count > 1)
            {
                lines[target_qubit_begin * 3 + 2                           ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[(target_qubit_begin + target_qubit_count - 1) * 3    ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 1].append("──┤").append(gate_strlen, ' ').append("├──");
                lines_length[target_qubit_begin * 3 + 2 ] += offset;
                lines_length[(target_qubit_begin + target_qubit_count - 1) * 3   ] += offset;
                lines_length[(target_qubit_begin + target_qubit_count - 1) * 3 + 1] += offset;
            }
            lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 2].append("  └").append(repeat(gate_strlen, "─")).append("┘  ");
            lines_length[(target_qubit_begin + target_qubit_count - 1) * 3 + 2] += offset;

            for (int i = target_qubit_begin + 1; i < target_qubit_begin + target_qubit_count - 1; ++i)
            {
                lines[i * 3    ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[i * 3 + 1].append("──┤").append(gate_strlen, ' ').append("├──");
                lines[i * 3 + 2].append("  │").append(gate_strlen, ' ').append("│  ");
                lines_length[i * 3    ] += offset;
                lines_length[i * 3 + 1] += offset;
                lines_length[i * 3 + 2] += offset;
            }
        }
        //Manage adding control qubits
        else
        {
            //Pad the complete circuit
            fill_circuit();

            lines[target_qubit_begin * 3 + 1].append("──┤").append(gate_str).append("├──");
            lines[control_qubit_begin * 3 + 1].append(repeat(gate_strmid + 3, "─")).append("█").append(repeat(gate_strlen - gate_strmid + 2, "─"));
            if (control_qubit_count > 1)
                lines[(control_qubit_begin + control_qubit_count - 1) * 3 + 1].append(
                    repeat(gate_strmid + 3, "─")).append("█").append(repeat(gate_strlen - gate_strmid + 2, "─"));

            lines[target_qubit_begin * 3    ].append("  ┌");
            lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 2].append("  └");
            lines[control_qubit_begin * 3    ].append(gate_strmid + 3, ' ');
            lines[control_qubit_begin * 3 + 2].append(gate_strmid + 3, ' ');
            if (control_qubit_count > 1)
            {
                lines[(control_qubit_begin + control_qubit_count - 1) * 3    ].append(gate_strmid + 3, ' ');
                lines[(control_qubit_begin + control_qubit_count - 1) * 3 + 2].append(gate_strmid + 3, ' ');
            }

            if (control_qubit_begin < target_qubit_begin)
            {
                lines[control_qubit_begin * 3 + 2].append("│");
                lines[target_qubit_begin * 3].append(repeat(gate_strmid, "─")).append("┴").append(
                    repeat(gate_strlen - gate_strmid - 1, "─"));
                lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 2].append(repeat(gate_strlen, "─"));
                if (control_qubit_count > 1)
                {
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3    ].append("│");
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3 + 2].append("│");
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3    ].append(gate_strlen - gate_strmid + 2, ' ');
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3 + 2].append(gate_strlen - gate_strmid + 2, ' ');
                }
                lines[control_qubit_begin * 3    ].append(gate_strlen - gate_strmid + 3, ' ');
                lines[control_qubit_begin * 3 + 2].append(gate_strlen - gate_strmid + 2, ' ');
            }
            else
            {
                lines[target_qubit_begin * 3].append(repeat(gate_strlen, "─"));
                lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 2].append(repeat(gate_strmid, "─")).append("┬").append(
                    repeat(gate_strlen - gate_strmid - 1, "─"));
                lines[control_qubit_begin * 3    ].append("│");

                if (control_qubit_count > 1)
                {
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3    ].append("│");
                    lines[control_qubit_begin * 3 + 2].append("│");
                    lines[(control_qubit_begin + control_qubit_count - 1) * 3    ].append(gate_strlen - gate_strmid + 2, ' ');
                    lines[control_qubit_begin * 3 + 2].append(gate_strlen - gate_strmid + 2, ' ');
                }
                lines[(control_qubit_begin + control_qubit_count - 1) * 3 + 2].append(gate_strlen - gate_strmid + 3, ' ');
                lines[control_qubit_begin * 3    ].append(gate_strlen - gate_strmid + 2, ' ');
            }
            
            lines[target_qubit_begin * 3    ].append("┐  ");
            lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 2].append("┘  ");

            lines_length[target_qubit_begin * 3    ] += offset;
            lines_length[target_qubit_begin * 3 + 1] += offset;
            lines_length[(target_qubit_begin + target_qubit_count - 1) * 3 + 2] += offset;
    
            lines_length[control_qubit_begin * 3    ] += offset;
            lines_length[control_qubit_begin * 3 + 1] += offset;
            lines_length[control_qubit_begin * 3 + 2] += offset;

            if (control_qubit_count > 1)
            {
                lines_length[(control_qubit_begin + control_qubit_count - 1) * 3    ] += offset;
                lines_length[(control_qubit_begin + control_qubit_count - 1) * 3 + 1] += offset;
                lines_length[(control_qubit_begin + control_qubit_count - 1) * 3 + 2] += offset;
            }

            for (int i = control_qubit_begin + 1; i < control_qubit_begin + control_qubit_count - 1; ++i)
            {
                lines[i * 3    ].append(gate_strmid + 3, ' ').append("│").append(gate_strlen - gate_strmid + 2, ' ');
                lines[i * 3 + 1].append(repeat(gate_strmid + 3, "─")).append("█").append(repeat(gate_strlen - gate_strmid + 2, "─"));
                lines[i * 3 + 2].append(gate_strmid + 3, ' ').append("│").append(gate_strlen - gate_strmid + 2, ' ');
                lines_length[i * 3    ] += offset;
                lines_length[i * 3 + 1] += offset;
                lines_length[i * 3 + 2] += offset;
            }

            if (target_qubit_count > 1)
            {
                lines[target_qubit_begin * 3 + 2                           ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[(target_qubit_begin + target_qubit_count - 1) * 3    ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[(target_qubit_begin + target_qubit_count - 1) * 3 + 1].append("──┤").append(gate_strlen, ' ').append("├──");
                lines_length[target_qubit_begin * 3 + 2 ] += offset;
                lines_length[(target_qubit_begin + target_qubit_count - 1) * 3   ] += offset;
                lines_length[(target_qubit_begin + target_qubit_count - 1) * 3 + 1] += offset;
            }

            for (int i = target_qubit_begin + 1; i < target_qubit_begin + target_qubit_count - 1; ++i)
            {
                lines[i * 3    ].append("  │").append(gate_strlen, ' ').append("│  ");
                lines[i * 3 + 1].append("──┤").append(gate_strlen, ' ').append("├──");
                lines[i * 3 + 2].append("  │").append(gate_strlen, ' ').append("│  ");
                lines_length[i * 3    ] += offset;
                lines_length[i * 3 + 1] += offset;
                lines_length[i * 3 + 2] += offset;
            }

            // Fill sections in between the control qubits and gates
            int begin_inter = control_qubit_begin < target_qubit_begin ? control_qubit_begin + control_qubit_count : target_qubit_begin + target_qubit_count;
            int end_inter = control_qubit_begin < target_qubit_begin ? target_qubit_begin : control_qubit_begin;
            for (int i = begin_inter; i < end_inter; ++i)
            {
                lines[i * 3    ].append(gate_strmid + 3, ' ').append("│").append(gate_strlen - gate_strmid + 2, ' ');
                lines[i * 3 + 1].append(repeat(gate_strmid + 3, "─")).append("┼").append(repeat(gate_strlen - gate_strmid + 2, "─"));
                lines[i * 3 + 2].append(gate_strmid + 3, ' ').append("│").append(gate_strlen - gate_strmid + 2, ' ');
                lines_length[i * 3    ] += offset;
                lines_length[i * 3 + 1] += offset;
                lines_length[i * 3 + 2] += offset;
            }

            // Align gates outside the range
            int min = control_qubit_begin < target_qubit_begin ? control_qubit_begin : target_qubit_begin;
            int max = (control_qubit_begin + control_qubit_count) < (target_qubit_begin + target_qubit_count) ?
                        target_qubit_begin + target_qubit_count : control_qubit_begin + control_qubit_count;
            int mid = (gate_strlen - padding.length() * 2) / 2;
            for (int i = 0; i < min; ++i)
            {
                lines[i * 3    ].append(mid, ' ');
                lines[i * 3 + 1].append(repeat(mid, "─"));
                lines[i * 3 + 2].append(mid, ' ');
                lines_length[i * 3    ] += mid;
                lines_length[i * 3 + 1] += mid;
                lines_length[i * 3 + 2] += mid;
            }
            for (int i = max; i < qubit_count; ++i)
            {
                lines[i * 3    ].append(mid, ' ');
                lines[i * 3 + 1].append(repeat(mid, "─"));
                lines[i * 3 + 2].append(mid, ' ');
                lines_length[i * 3    ] += mid;
                lines_length[i * 3 + 1] += mid;
                lines_length[i * 3 + 2] += mid;
            }
        }
       
        if (lines_length[target_qubit_begin * 3] > current_circuit_length)
            current_circuit_length = lines_length[target_qubit_begin * 3];

        gates_current_begin = gates_current_end + 1;
    }

    //Fill circuit till end
    fill_circuit();

    //Add newlines
    auto total_length = 0;
    for (auto& line : lines)
    {
        line.append("\n");
        total_length += line.length();
    }

    std::string out;
    out.reserve(total_length + 1);
    out.append(1, '\n');

    //Append all lines into output string
    for (const auto& line : lines)
        out.append(line);

    return out;
}

std::string gen_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator)
{
    if (circuit.qubit_count() != simulator.qubit_count())
        throw std::invalid_argument{"Circuit and simulator qubit count must match"};

    //Add number of qubits
    const int qubits = circuit.qubit_count();

    std::string circuit_representation = std::to_string(qubits) + "\n";

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

        circuit_representation.append(std::to_string(i) + "," + std::to_string(val) + ";");
    }

    //Add circuit gates
    circuit_representation.append("\n");
    circuit_representation.append(circuit.representation());

    return parse_circuit_representation(circuit_representation);
}

void print_circuit_text_image(const QCircuit& circuit, const QSimulator& simulator)
{
    std::cout << gen_circuit_text_image(circuit, simulator) << std::endl;
}
}