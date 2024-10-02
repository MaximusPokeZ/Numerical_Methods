cmake_minimum_required(VERSION 3.27)

add_subdirectory(test)

add_library(Rotate src/rotate_method.cpp)

target_include_directories(
        Rotate
        PUBLIC
        ./include
)

target_link_libraries(
        Rotate
        PUBLIC
        Matrix
)