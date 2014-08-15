cmake_minimum_required (VERSION 2.8)
project(verified-math)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory("${PROJECT_SOURCE_DIR}/gtest-1.7.0")
include_directories("${PROJECT_SOURCE_DIR}/gtest-1.7.0/include")

add_subdirectory("${PROJECT_SOURCE_DIR}/checkpp")
include_directories("${PROJECT_SOURCE_DIR}/checkpp/include")

add_library(verified_math SHARED
 src/main/vec3.cpp
)

add_executable(test_vec3
  src/test/test_vec3.cpp
)
target_link_libraries(test_vec3 gtest_main checkpp)
