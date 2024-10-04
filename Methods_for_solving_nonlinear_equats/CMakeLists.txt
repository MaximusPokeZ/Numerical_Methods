cmake_minimum_required(VERSION 3.27)

add_subdirectory(test)

add_library(Methods src/methods.cpp)

target_include_directories(
        Methods
        PUBLIC
        ./include
)


