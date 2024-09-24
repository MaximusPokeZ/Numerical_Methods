cmake_minimum_required(VERSION 3.22)

project(ThomasMethod)

add_subdirectory(test)
add_library(ThomasMethod src/method.cpp)

target_include_directories(
        ThomasMethod
        PUBLIC
        ./include)

target_link_libraries(
        ThomasMethod
        PUBLIC
        Matrix)

set_target_properties(
        ThomasMethod PROPERTIES
        LANGUAGES CXX
        LINKER_LANGUAGE CXX
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Thomas method")
