# Copyright (c) 2022, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(afQuantumSimulator
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)

find_package(ArrayFire REQUIRED)
find_package(nlopt REQUIRED)

add_library(afquantum
                src/quantum.cpp
                src/quantum_algo.cpp
                src/quantum_gates.cpp
                src/quantum_visuals.cpp
                src/utils.cpp
                include/quantum.h
                include/quantum_algo.h
                include/quantum_gates.h
                include/quantum_visuals.h
                include/utils.h)

target_include_directories(afquantum PUBLIC include)

set(AF_BUILD_BACKEND ArrayFire::af CACHE STRING "Arrayfire backend type")
target_link_libraries(afquantum PRIVATE NLopt::nlopt PUBLIC ${AF_BUILD_BACKEND})

option(AQS_BUILD_EXAMPLES   "Build AQS Examples" ON)
option(AQS_BUILD_TESTS      "Build AQS Tests"    ON)
option(AQS_BUILD_BENCHMARKS "Build benchmarks"   OFF)

if (${AQS_BUILD_EXAMPLES})
        add_subdirectory(examples)
endif()
if (${AQS_BUILD_TESTS})
        add_subdirectory(test)
endif()
if (${AQS_BUILD_BENCHMARKS})
        add_subdirectory(benchmark)
endif()