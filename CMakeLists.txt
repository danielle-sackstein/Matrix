cmake_minimum_required(VERSION 3.8)
project(ex3)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp Matrix.cpp Matrix.hpp Exceptions.cpp Exceptions.h Complex.h
        Complex.cpp GenericMatrixDriver.cpp ParallelChecker.cpp TimeChecker.cpp)
include_directories(".")
add_executable(ex3 ${SOURCE_FILES})