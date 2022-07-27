# Copyright (c) 2022, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(afQuantumSim-Examples
        VERSION 1.0
        LANGUAGES CXX)

add_executable(example examples.cpp)
target_link_libraries(example PRIVATE afquantum)

add_executable(helloworld helloworld.cpp)
target_link_libraries(helloworld PRIVATE afquantum)

add_executable(superposition superposition.cpp)
target_link_libraries(superposition PRIVATE afquantum)

add_executable(entanglement entanglement.cpp)
target_link_libraries(entanglement PRIVATE afquantum)

add_executable(basis_change basis_change.cpp)
target_link_libraries(basis_change PRIVATE afquantum)

add_executable(classical_gates classical_gates.cpp)
target_link_libraries(classical_gates PRIVATE afquantum)

add_executable(classic_2bit_adder classic_2bit_adder.cpp)
target_link_libraries(classic_2bit_adder PRIVATE afquantum)

add_executable(fourier_transform fourier_transform.cpp)
target_link_libraries(fourier_transform PRIVATE afquantum)

add_executable(qft_adder qft_adder.cpp)
target_link_libraries(qft_adder PRIVATE afquantum)

add_executable(phase_estimation phase_estimation.cpp)
target_link_libraries(phase_estimation PRIVATE afquantum)

add_executable(grover_search grover_search.cpp)
target_link_libraries(grover_search PRIVATE afquantum)

add_executable(quantum_counting quantum_counting.cpp)
target_link_libraries(quantum_counting PRIVATE afquantum)

add_executable(quantum_teleportation quantum_teleportation.cpp)
target_link_libraries(quantum_teleportation PRIVATE afquantum)

add_executable(vqe vqe.cpp)
target_link_libraries(vqe PRIVATE afquantum)

add_executable(vqr vqr.cpp)
target_link_libraries(vqr PRIVATE afquantum NLopt::nlopt)