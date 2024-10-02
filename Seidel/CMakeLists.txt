cmake_minimum_required(VERSION 3.27)

project(seidel)

add_subdirectory(test)

add_library(Seidel src/functions.cpp)

target_include_directories(
        Seidel
        PUBLIC
        ./include
)

target_link_libraries(
        Seidel
        PUBLIC
        Matrix
)

set_target_properties(
        Seidel PROPERTIES
        LANGUAGES CXX
        LINKER_LANGUAGE CXX
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Seidel method"
)