cmake_minimum_required(VERSION 3.6)
project(Particle-Sim-MPI)
File(GLOB SRC_FILES src/*.cpp)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=gnu++11 -lm -o3")

message(src: ${SRC_FILES})
add_executable(Particle-Sim-MPI ${SRC_FILES} )

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(Particle-Sim-MPI ${MPI_LIBRARIES})