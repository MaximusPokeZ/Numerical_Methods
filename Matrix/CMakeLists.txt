cmake_minimum_required(VERSION 3.22)

project(Matrix)

add_library(Matrix src/Matrix.cpp)

target_include_directories(Matrix PUBLIC ./include)

set_target_properties(
        Matrix PROPERTIES
        LANGUAGES CXX
        LINKER_LANGUAGE CXX
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Matrix interface library"
)
