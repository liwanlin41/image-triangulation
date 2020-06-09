cmake_minimum_required(VERSION 3.10)
project(CoarseTriangulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

add_executable(CoarseTriangulation main.cpp)

target_link_libraries(CoarseTriangulation PUBLIC src)

target_include_directories(CoarseTriangulation PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")
