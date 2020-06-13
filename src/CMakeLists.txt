# from https://stackoverflow.com/questions/3201154
file(GLOB src_files "*.cpp")
#add_library(src point.cpp segment.cpp matrix.cpp triangle.cpp)
add_library(src ${src_files})

find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})
link_directories(${X11_LIBRARIES})

target_include_directories(src INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
