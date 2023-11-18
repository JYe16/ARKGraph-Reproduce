cmake_minimum_required(VERSION 3.26)
project(ARKGraph_Reproduce)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(ARKGraph_Reproduce main.cpp)
