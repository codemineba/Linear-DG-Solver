cmake_minimum_required(VERSION 3.26)
project(DG)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -O3")

# 添加外部库
add_subdirectory(${CMAKE_SOURCE_DIR}/external/SciCalUtils)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/mesh_structure)

# 添加子模块
add_subdirectory(src/LinearDGSolver2D)
add_subdirectory(src/LinearDGSolver3D)
add_subdirectory(demo)
