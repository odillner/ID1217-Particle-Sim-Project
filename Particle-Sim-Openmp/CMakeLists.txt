cmake_minimum_required(VERSION 3.6)
project(Particle-Sim-Openmp)

File(GLOB SRC_FILES src/*.cpp)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=gnu++11 -lm -o3 -fopenmp" )

message(src: ${SRC_FILES})
add_executable(Particle-Sim-Openmp ${SRC_FILES})