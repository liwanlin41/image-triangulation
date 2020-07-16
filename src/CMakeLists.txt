# from https://stackoverflow.com/questions/3201154
file(GLOB src_files "*.cpp" "*.cu")
add_library(src ${src_files})

#add_subdirectory("../deps/polyscope" polyscope)

#include_directories("../deps/polyscope" polyscope)
#target_include_directories(src PUBLIC "../deps/polyscope/include")
#target_include_directories(src PUBLIC polyscope)

# get matlab
find_package(Matlab REQUIRED)
include_directories("${Matlab_ROOT_DIR}/extern/include")

target_link_libraries(src PUBLIC ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so)
target_link_libraries(src PUBLIC ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so)


#target_include_directories(src PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/polyscope")
target_include_directories(src INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
#target_link_libraries(src polyscope)
