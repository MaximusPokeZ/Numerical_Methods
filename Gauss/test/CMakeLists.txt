cmake_minimum_required(VERSION 3.22)

project(gauss_method_test)

add_executable(gauss_method_test main.cpp)

target_link_libraries(
        gauss_method_test
        PUBLIC
        Matrix
)

target_link_libraries(
        gauss_method_test
        PUBLIC
        Gauss
)

