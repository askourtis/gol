cmake_minimum_required(VERSION 3.15...3.30)

project(
  libgol
  VERSION 1.0
  LANGUAGES C CXX CUDA
)


file(GLOB SOURCES "${CMAKE_CURRENT_LIST_DIR}/*.c"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cpp"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cc"
                  "${CMAKE_CURRENT_LIST_DIR}/*.cu"
)
message("${PROJECT_NAME}")
add_library(
  "${PROJECT_NAME}"
  "${SOURCES}"
)

target_include_directories("${PROJECT_NAME}" INTERFACE "${CMAKE_CURRENT_LIST_DIR}/../inc")

