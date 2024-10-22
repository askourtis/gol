cmake_minimum_required(VERSION 3.15...3.30)

project(
  testgol
  VERSION 1.0
  LANGUAGES C CXX CUDA
)


file(GLOB SOURCES "${CMAKE_CURRENT_LIST_DIR}/*.c"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cpp"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cc"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cu"
)

add_executable(
  "${PROJECT_NAME}"
  "${SOURCES}"
)


#add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../src")

target_link_libraries("${PROJECT_NAME}" libgol)

