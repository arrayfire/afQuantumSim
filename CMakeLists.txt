# Copyright (c) 2022, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(afQuantumSimulator)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)

find_package(ArrayFire REQUIRED)
find_package(OpenCL QUIET)
find_package(CUDA QUIET)
find_package(nlopt REQUIRED)

if(${ArrayFire_CPU_FOUND})
  add_executable(test_scientific_cpu tests.cpp quantum.cpp utils.cpp quantum_visuals.cpp quantum_gates.cpp)
  add_executable(example_scientific_cpu examples.cpp quantum.cpp utils.cpp quantum_algo.cpp quantum_gates.cpp quantum_visuals.cpp)
  add_executable(stress_scientific_cpu quantum.cpp stress_test.cpp utils.cpp)
  add_executable(simulation_cpu quantum.cpp utils.cpp simulation.cpp quantum_visuals.cpp)
  target_link_libraries(stress_scientific_cpu ArrayFire::afcpu)
  target_link_libraries(example_scientific_cpu ArrayFire::afcpu)
  target_link_libraries(test_scientific_cpu ArrayFire::afcpu)
  target_link_libraries(simulation_cpu ArrayFire::afcpu NLopt::nlopt)
endif()

if(${ArrayFire_CUDA_FOUND})
  add_executable(test_scientific_cuda tests.cpp quantum.cpp utils.cpp quantum_visuals.cpp quantum_gates.cpp)
  add_executable(example_scientific_cuda examples.cpp quantum.cpp utils.cpp quantum_algo.cpp quantum_gates.cpp quantum_visuals.cpp)
  add_executable(stress_scientific_cuda quantum.cpp stress_test.cpp utils.cpp)
  add_executable(simulation_cuda quantum.cpp utils.cpp simulation.cpp quantum_visuals.cpp)
  target_link_libraries(stress_scientific_cuda ArrayFire::afcuda)
  target_link_libraries(example_scientific_cuda ArrayFire::afcuda)
  target_link_libraries(test_scientific_cuda ArrayFire::afcuda)
  target_link_libraries(simulation_cuda ArrayFire::afcuda NLopt::nlopt)
endif()

if(${ArrayFire_OpenCL_FOUND})
  add_executable(test_scientific_opencl tests.cpp quantum.cpp utils.cpp quantum_visuals.cpp quantum_gates.cpp)
  add_executable(example_scientific_opencl examples.cpp quantum.cpp utils.cpp quantum_algo.cpp quantum_gates.cpp quantum_visuals.cpp)
  add_executable(stress_scientific_opencl quantum.cpp stress_test.cpp utils.cpp)
  add_executable(simulation_opencl quantum.cpp utils.cpp simulation.cpp quantum_visuals.cpp)
  target_link_libraries(stress_scientific_opencl ArrayFire::afopencl)
  target_link_libraries(example_scientific_opencl ArrayFire::afopencl)
  target_link_libraries(test_scientific_opencl ArrayFire::afopencl)
  target_link_libraries(simulation_opencl ArrayFire::afopencl NLopt::nlopt)
endif()