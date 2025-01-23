cmake_minimum_required(VERSION 3.10)
project(QP_based_motion_filter)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# qpOASES
add_subdirectory(qpOASES)
include_directories(qpOASES/include)

# QP based motion filter
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

add_executable(QP_based_motion_filter  ${SOURCES})

target_link_libraries(QP_based_motion_filter qpOASES)
