cmake_minimum_required(VERSION 3.10)
project(crop)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(crop crop.cpp)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})
find_package(X11 REQUIRED)
include_directories(${X11_LIBRARIES})
link_directories(${X11_LIBRARIES})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(crop PUBLIC Threads::Threads)

target_link_libraries(crop PUBLIC ${PNG_LIBRARY})
target_link_libraries(crop PUBLIC ${JPEG_LIBRARY})
target_link_libraries(crop PUBLIC ${X11_LIBRARIES})
