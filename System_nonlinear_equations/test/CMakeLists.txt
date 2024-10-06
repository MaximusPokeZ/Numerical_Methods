cmake_minimum_required(VERSION 3.27)

add_executable(system_nonlinear_test main.cpp)

target_link_libraries(
        system_nonlinear_test
        PUBLIC
        System_nonlinear
)

target_link_libraries(
        system_nonlinear_test
        PUBLIC
        Matrix
)
