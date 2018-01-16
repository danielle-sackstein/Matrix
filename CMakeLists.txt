cmake_minimum_required(VERSION 3.8)
project(ex3)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp matrix.cpp matrix.hpp MatrixException.cpp MatrixException.h)
add_executable(ex3 ${SOURCE_FILES})