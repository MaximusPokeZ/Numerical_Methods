cmake_minimum_required(VERSION 3.22)
project(Numerical_Methods VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(Matrix)
add_subdirectory(Thomas_method)
add_subdirectory(Gauss)
add_subdirectory(Seidel)
add_subdirectory(Rotation_method)