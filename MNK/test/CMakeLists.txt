add_executable(mnk main.cpp)

target_link_directories(mnk PRIVATE ../include)

target_link_libraries(
        mnk
        PUBLIC
        Matrix
)

target_link_libraries(
        mnk
        PUBLIC
        Gauss
)