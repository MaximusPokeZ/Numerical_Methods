cmake_minimum_required(VERSION 3.27)

project(seidel_method_test)

add_executable(simple_iter simple_iter.cpp)
add_executable(seidel seidel.cpp)

target_link_libraries(
        simple_iter
        PUBLIC
        Matrix
)

target_link_libraries(
        simple_iter
        PUBLIC
        Seidel
)

target_link_libraries(
        seidel
        PUBLIC
        Matrix
)

target_link_libraries(
        seidel
        PUBLIC
        Seidel
)