cmake_minimum_required(VERSION 3.27)

project(rotate_method_test)

add_executable(rotate_method_test main.cpp)

target_link_libraries(
        rotate_method_test
        PUBLIC
        Matrix
)

target_link_libraries(
        rotate_method_test
        PUBLIC
        Rotate
)