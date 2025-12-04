# add cmake version, project name and standard
cmake_minimum_required(VERSION 3.20)
project(mechsystem)

set (CMAKE_CXX_STANDARD 20)

# include the head file in the root 
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# file need to be compiled
add_executable (test_mass_spring mass_spring.cpp ${CMAKE_SOURCE_DIR}/src)



find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)


#use pybind to create python module and privately:will not influence other file
pybind11_add_module(mass_spring
    bind_mass_spring.cpp
)

target_include_directories(mass_spring PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)


