cmake_minimum_required(VERSION 3.22)

project(thomas_method_test)

add_executable(thomas_method_test main.cpp)

target_link_libraries(
        thomas_method_test
        PUBLIC
        Matrix
)

target_link_libraries(
        thomas_method_test
        PUBLIC
        ThomasMethod
)

set_target_properties(
        thomas_method_test PROPERTIES
        LANGUAGES CXX
        LINKER_LANGUAGE CXX
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Thomas method test")