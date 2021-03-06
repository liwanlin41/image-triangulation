# taken from polyscope CMakeLists

cmake_minimum_required(VERSION 3.10)

# Set up testing
enable_testing()

project(triangulation-tests LANGUAGES CUDA CXX)

find_package(CUDA 10.1 REQUIRED)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Configure the compiler

#if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#
#  # using Clang (linux or apple) or GCC
#  message("Using clang/gcc compiler flags")
#  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -g3")
#  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
#  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")
#
#  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#    message("Setting clang-specific options")
#    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
#    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
#  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    message("Setting gcc-specific options")
#    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
#    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
#  endif()
#  
#  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
#  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
#  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
#
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#  # using Visual Studio C++
#  message("Using Visual Studio compiler flags")
#  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
#  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
#  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
#  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
#  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
#  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
#  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
#  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
#
#  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
#  add_definitions (-DNOMINMAX)
#else()
#  # unrecognized
#  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
#endif()



### Download and unpack googletest at configure time
# (from the googletest docs)
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

set_property(TARGET gtest PROPERTY CXX_STANDARD 17)

# include X11 for CIMg
find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})
link_directories(${X11_LIBRARIES})

# Build the tests
set(TEST_SRCS structureTest.cpp integratorTest.cpp gradientTest.cpp)

add_executable(triangulation-test "${TEST_SRCS}")
target_include_directories(triangulation-test PUBLIC "${gtest_SOURCE_DIR}/include")
target_link_libraries(triangulation-test gtest_main src ${X11_LIBRARIES})
