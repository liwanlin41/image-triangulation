cmake_minimum_required(VERSION 3.10)
project(CoarseTriangulation LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

set_target_properties(src PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(CoarseTriangulation main.cpp)
#set_property(TARGET CoarseTriangulation PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# from https://stackoverflow.com/questions/38903996/
# include png libraries
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

# include jpeg libraries
find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

# include polyscope
add_subdirectory(deps/polyscope)

# include matlab
find_package(Matlab REQUIRED)
include_directories("${Matlab_ROOT_DIR}/extern/include")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(CoarseTriangulation PUBLIC Threads::Threads)
target_link_libraries(CoarseTriangulation PUBLIC ${PNG_LIBRARY})
target_link_libraries(CoarseTriangulation PUBLIC ${JPEG_LIBRARY})
target_link_libraries(CoarseTriangulation PUBLIC src)
target_link_libraries(CoarseTriangulation PUBLIC polyscope)
target_link_libraries(CoarseTriangulation PUBLIC ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so)
target_link_libraries(CoarseTriangulation PUBLIC ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so)

target_include_directories(CoarseTriangulation PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")
