cmake_minimum_required(VERSION 3.15...3.30)

project(GameOfLifeProject)



set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src")
set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

add_subdirectory(${SRC_DIR})
add_subdirectory(${TEST_DIR})

