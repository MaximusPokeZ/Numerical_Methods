add_executable(spline main.cpp)

target_link_directories(spline PRIVATE ../include)

target_link_libraries(
        spline
        PUBLIC
        Matrix
)

target_link_libraries(
        spline
        PUBLIC
        Gauss
)