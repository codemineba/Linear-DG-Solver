cmake_minimum_required(VERSION 3.26)
project(DG)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -O3")

# include(FetchContent)
# # 下载并添加 SciCalUtils
# FetchContent_Declare(
#     SciCalUtils
#     GIT_REPOSITORY https://github.com/codemineba/SciCalUtils.git
#     GIT_TAG master
#     SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/SciCalUtils
# )
# FetchContent_MakeAvailable(SciCalUtils)  # 自动执行 add_subdirectory

# # 下载并添加 mesh_structure
# FetchContent_Declare(
#     mesh_structure
#     GIT_REPOSITORY https://github.com/codemineba/mesh_structure.git
#     GIT_TAG master
#     SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/mesh_structure
# )
# FetchContent_MakeAvailable(mesh_structure)

# 下载并添加 Pybind11
# FetchContent_Declare(
#     pybind11
#     GIT_REPOSITORY https://github.com/pybind/pybind11.git
#     GIT_TAG master
#     SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/pybind11
# )
# FetchContent_MakeAvailable(pybind11)  # 自动执行 add_subdirectory


add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/SciCalUtils)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/mesh_structure)

# 添加子模块
add_subdirectory(src/LinearDGSolver2D)
add_subdirectory(src/LinearDGSolver3D)
add_subdirectory(pymodule)
add_subdirectory(demo)
