cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
find_package(qpp REQUIRED)

project(AQS-Benchmark)
add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark afquantum)

add_executable(qpp_benchmark qpp_benchmark.cpp)
target_link_libraries(qpp_benchmark PUBLIC ${QPP_LINK_DEPS} libqpp)