cmake_minimum_required(VERSION 2.8.3)
find_package( OpenCV REQUIRED )
include_directories(include)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

add_library(network src/network.cpp)
add_library(kmeans src/kmeans.cpp)
add_library(pcnn_simple src/pcnn_simple.cpp)
add_library(Line_detector_class src/Line_detector.cpp)

# add_executable(test src/main.cpp)
add_executable(line_segmentation src/line_segmentation.cpp)

target_link_libraries(pcnn_simple network)
target_link_libraries(line_segmentation pcnn_simple network Line_detector_class kmeans ${OpenCV_LIBS})