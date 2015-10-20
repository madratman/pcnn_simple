cmake_minimum_required(VERSION 2.8.3)
include_directories(include)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

add_library(network src/network.cpp)
add_library(pcnn src/pcnn.cpp)

add_executable(main src/main.cpp)

target_link_libraries(pcnn network)
target_link_libraries(main pcnn network)
