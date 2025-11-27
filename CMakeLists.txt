cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)


include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)

add_subdirectory (mechsystem)

add_executable (test_ode demos/test_ode.cpp)
target_link_libraries (test_ode PUBLIC nanoblas)

add_executable(ex18_5 Exercise2/ex18_5.cpp)
target_link_libraries(ex18_5 PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

add_executable(Explicit_Euler
    Exercise1/ex17.2.2/Explicit_Euler.cpp
)
target_link_libraries(Explicit_Euler PUBLIC nanoblas)

add_executable(Improved_Euler
    Exercise1/ex17.2.2/Improved_Euler.cpp
)
target_link_libraries(Improved_Euler PUBLIC nanoblas)


add_executable(Implicit_Euler Exercise1/ex17.4.1/Implicit_Euler.cpp)
target_link_libraries(Implicit_Euler PUBLIC nanoblas)

add_executable(Crank_Nicolson Exercise1/ex17.4.1/Crank_Nicolson.cpp)
target_link_libraries(Crank_Nicolson PUBLIC nanoblas)
