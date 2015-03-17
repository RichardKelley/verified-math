cmake_minimum_required (VERSION 2.8)
project(verified-math)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory("${PROJECT_SOURCE_DIR}/gtest-1.7.0")
include_directories("${PROJECT_SOURCE_DIR}/gtest-1.7.0/include")

add_subdirectory("${PROJECT_SOURCE_DIR}/checkpp")
include_directories("${PROJECT_SOURCE_DIR}/checkpp/include")

add_library(verified_math SHARED
 src/main/vec3.cpp
 src/main/vec4.cpp
)

add_executable(test_vec3
  src/test/test_vec3.cpp
)
target_link_libraries(test_vec3 gtest_main checkpp)

add_executable(test_vec4
  src/test/test_vec4.cpp
)
target_link_libraries(test_vec4 gtest_main checkpp)

add_executable(test_mat33
  src/test/test_mat33.cpp
)
target_link_libraries(test_mat33 gtest_main checkpp)

add_executable(test_mat44
  src/test/test_mat44.cpp
)
target_link_libraries(test_mat44 gtest_main checkpp)
