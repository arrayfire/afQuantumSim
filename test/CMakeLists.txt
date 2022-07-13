# Copyright (c) 2022, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(afQuantumSim-Tests
        VERSION 1.0
        LANGUAGES CXX)

if(${ArrayFire_CPU_FOUND})
    add_executable(test_scientific_cpu tests.cpp ${AQS_SOURCES})
    add_executable(stress_scientific_cpu stress_test.cpp stress_test.cpp ${AQS_SOURCES})
    target_link_libraries(stress_scientific_cpu ArrayFire::afcpu NLopt::nlopt)
    target_link_libraries(test_scientific_cpu ArrayFire::afcpu NLopt::nlopt)
    target_include_directories(stress_scientific_cpu PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(test_scientific_cpu PUBLIC ${AQS_HEADER_DIR})
endif()

if(${ArrayFire_CUDA_FOUND})
    add_executable(test_scientific_cuda tests.cpp ${AQS_SOURCES})
    add_executable(stress_scientific_cuda stress_test.cpp ${AQS_SOURCES})
    target_link_libraries(stress_scientific_cuda ArrayFire::afcuda NLopt::nlopt)
    target_link_libraries(test_scientific_cuda ArrayFire::afcuda NLopt::opt)
    target_include_directories(stress_scientific_cuda PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(test_scientific_cuda PUBLIC ${AQS_HEADER_DIR})
endif()

if(${ArrayFire_OpenCL_FOUND})
    add_executable(test_scientific_opencl tests.cpp ${AQS_SOURCES})
    add_executable(stress_scientific_opencl stress_test.cpp ${AQS_SOURCES})
    target_link_libraries(stress_scientific_opencl ArrayFire::afopencl NLopt::nlopt)
    target_link_libraries(test_scientific_opencl ArrayFire::afopencl NLopt::nlopt)
    target_include_directories(stress_scientific_opencl PUBLIC ${AQS_HEADER_DIR})
    target_include_directories(test_scientific_opencl PUBLIC ${AQS_HEADER_DIR})
endif()