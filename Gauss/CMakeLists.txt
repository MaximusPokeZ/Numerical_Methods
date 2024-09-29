cmake_minimum_required(VERSION 3.27)
project(Gauss)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(test)

add_library(Gauss src/gauss.cpp)

target_include_directories(
        Gauss
        PUBLIC
        ./include
)

target_link_libraries(
        Gauss
        PUBLIC
        Matrix
)

set_target_properties(
        Gauss PROPERTIES
        LANGUAGES CXX
        LINKER_LANGUAGE CXX
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Gauss method"
)


