cmake_minimum_required(VERSION 3.10)
project(mesh_structure)

# 创建 triangle_mesh 库
add_library(triangle_mesh src/TriangleMesh.cpp)
target_include_directories(triangle_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 创建 tetrahedron_mesh 库
add_library(tetrahedron_mesh src/TetrahedronMesh.cpp)
target_include_directories(tetrahedron_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
