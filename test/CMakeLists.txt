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

add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE afquantum)

add_executable(stress_test stress_test.cpp)
target_link_libraries(stress_test PRIVATE afquantum)