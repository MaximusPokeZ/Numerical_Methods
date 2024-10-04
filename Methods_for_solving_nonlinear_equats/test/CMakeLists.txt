cmake_minimum_required(VERSION 3.27)

project(nonlinear_methods_test)

add_executable(nonlinear_methods_test main.cpp)

target_link_libraries(
        nonlinear_methods_test
        PUBLIC
        Methods
)