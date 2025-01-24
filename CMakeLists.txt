cmake_minimum_required(VERSION 3.15...3.30)

project(GameOfLifeProject)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CXX_STANDARD 23)
set(CXX_STANDARD_REQUIRED On)

set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src")
set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/test")
set(BMARK_DIR "${CMAKE_CURRENT_LIST_DIR}/benchmark")

add_subdirectory(${SRC_DIR})
add_subdirectory(${TEST_DIR})
add_subdirectory(${BMARK_DIR})
