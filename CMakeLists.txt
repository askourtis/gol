cmake_minimum_required(VERSION 3.15...3.30)

project(GameOfLifeProject)



set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/src")

add_subdirectory(${SRC_DIR})

add_executable(
  testgol
  test/test.c
)

target_link_libraries(testgol PRIVATE libgol)

