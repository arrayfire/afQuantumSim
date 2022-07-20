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

if(${ArrayFire_CPU_FOUND})
    add_executable(example_scientific_cpu examples.cpp ${AQS_SOURCES})
    add_executable(vqe_example_cpu vqe_example.cpp ${AQS_SOURCES} ${AQS_HEADERS})
    add_executable(vqr_example_cpu vqr_example.cpp ${AQS_SOURCES})
    target_link_libraries(example_scientific_cpu ArrayFire::afcpu NLopt::nlopt)
    target_link_libraries(vqe_example_cpu ArrayFire::afcpu NLopt::nlopt)
    target_link_libraries(vqr_example_cpu ArrayFire::afcpu NLopt::nlopt)
    target_include_directories(example_scientific_cpu PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqe_example_cpu PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqr_example_cpu PUBLIC ${AQS_HEADER_DIR})
endif()

if(${ArrayFire_CUDA_FOUND})
    add_executable(example_scientific_cuda examples.cpp ${AQS_SOURCES})
    add_executable(vqe_example_cuda vqe_example.cpp ${AQS_SOURCES})
    add_executable(vqr_example_cuda vqr_example.cpp ${AQS_SOURCES})
    target_link_libraries(example_scientific_cuda ArrayFire::afcuda NLopt::nlopt)
    target_link_libraries(vqe_example_cuda ArrayFire::afcuda NLopt::nlopt)
    target_link_libraries(vqr_example_cuda ArrayFire::afcuda NLopt::nlopt)
    target_include_directories(example_scientific_cuda PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqe_example_cuda PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqr_example_cuda PUBLIC ${AQS_HEADER_DIR})
endif()

if(${ArrayFire_OpenCL_FOUND})
    add_executable(example_scientific_opencl examples.cpp ${AQS_SOURCES})
    add_executable(vqe_example_opencl vqe_example.cpp ${AQS_SOURCES})
    add_executable(vqr_example_opencl vqr_example.cpp ${AQS_SOURCES})
    target_link_libraries(example_scientific_opencl ArrayFire::afopencl NLopt::nlopt)
    target_link_libraries(vqe_example_opencl ArrayFire::afopencl NLopt::nlopt)
    target_link_libraries(vqr_example_opencl ArrayFire::afopencl NLopt::nlopt)
    target_include_directories(example_scientific_opencl PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqe_example_opencl PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(vqr_example_opencl PUBLIC ${AQS_HEADER_DIR})
endif()