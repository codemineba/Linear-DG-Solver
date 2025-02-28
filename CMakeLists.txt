cmake_minimum_required(VERSION 3.26)
project(2D_LE)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -O3")


# 添加子模块
add_subdirectory(src/LinearDGSolver2D)
add_subdirectory(src/LinearDGSolver3D)
add_subdirectory(src/mesh_structure)
add_subdirectory(src/Utils)
add_subdirectory(demo)
