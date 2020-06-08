add_library(src point.cpp segment.cpp matrix.cpp triangle.cpp)

target_include_directories(src INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
