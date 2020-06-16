cmake_minimum_required(VERSION 3.10)
project(CoarseTriangulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

add_executable(CoarseTriangulation main.cpp)

# from https://stackoverflow.com/questions/38903996/
# include png libraries
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

# include jpeg libraries
find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

# include polyscope
add_subdirectory(deps/polyscope)

#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)
#target_link_libraries(CoarseTriangulation PUBLIC Threads::Threads)
#target_link_libraries(CoarseTriangulation PUBLIC X11)
target_link_libraries(CoarseTriangulation PUBLIC ${PNG_LIBRARY})
target_link_libraries(CoarseTriangulation PUBLIC ${JPEG_LIBRARY})
target_link_libraries(CoarseTriangulation PUBLIC src)
target_link_libraries(CoarseTriangulation PUBLIC polyscope)

target_include_directories(CoarseTriangulation PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")
