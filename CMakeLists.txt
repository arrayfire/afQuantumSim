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
find_package(OpenCL QUIET)
find_package(CUDA QUIET)
find_package(nlopt REQUIRED)

add_subdirectory(examples)
add_subdirectory(test)

set(AQS_PROJECT_DIR "${PROJECT_SOURCE_DIR}")
message(${AQS_PROJECT_DIR})
set(AQS_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/quantum.cpp"
                "${CMAKE_CURRENT_LIST_DIR}/src/quantum_algo.cpp"
                "${CMAKE_CURRENT_LIST_DIR}/src/quantum_gates.cpp"
                "${CMAKE_CURRENT_LIST_DIR}/src/quantum_visuals.cpp"
                "${CMAKE_CURRENT_LIST_DIR}/src/utils.cpp"
                CACHE INTERNAL "")

set(AQS_HEADER_DIR "${CMAKE_CURRENT_LIST_DIR}/include" CACHE INTERNAL "")

set(AQS_HEADERS "${CMAKE_CURRENT_LIST_DIR}/include/quantum.h"
                "${CMAKE_CURRENT_LIST_DIR}/include/quantum_algo.h"
                "${CMAKE_CURRENT_LIST_DIR}/include/quantum_gates.h"
                "${CMAKE_CURRENT_LIST_DIR}/include/quantum_visuals.h"
                "${CMAKE_CURRENT_LIST_DIR}/include/utils.h"
                CACHE INTERNAL "")