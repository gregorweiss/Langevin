cmake_minimum_required(VERSION 3.15)
project(Langevin)

set(CMAKE_CXX_STANDARD 14)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -D_GLIBCXX_PARALLEL")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized=vec.info")
# Eigen Linear Algebra
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/eigen3 -L/usr/lib/x86_64-linux-gnu/eigen3 ")

add_executable(Langevin main.cpp
        src/ArgParse.cpp src/ArgParse.h
        src/Argument.cpp src/Argument.h
        src/common.h
        src/Force.cpp src/Force.h
        src/HIcouple.cpp src/HIcouple.h
        src/Particle.cpp src/Particle.h
        src/System.cpp src/System.h)
