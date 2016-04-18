cmake_minimum_required(VERSION 3.3)
project(MetricSpaceExample)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES main.cpp MetricTree.hpp FaSST.hpp FaSSTGating.hpp)
add_executable(MetricSpaceExample ${SOURCE_FILES})