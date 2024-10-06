cmake_minimum_required(VERSION 3.27)

add_subdirectory(test)

add_library(System_nonlinear ./src/method.cpp)

target_include_directories(
        System_nonlinear
        PUBLIC
        ./include
)

target_link_libraries(
        System_nonlinear
        PUBLIC
        Matrix
)
