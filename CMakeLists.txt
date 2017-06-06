cmake_minimum_required(VERSION 3.6)
project(discretelogarithms)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c)
add_executable(discretelogarithms ${SOURCE_FILES})
add_compile_options(O2)

target_link_libraries(discretelogarithms m)
target_link_libraries(discretelogarithms gmp)
