include_directories(${CMAKE_SOURCE_DIR}/include/mesh_structure)
include_directories(${CMAKE_SOURCE_DIR}/include/Utils)

# triangle_mesh_structure
add_library(triangle_mesh
    TriangleMesh.cpp
)

# tetrahedron_mesh_structure
add_library(tetrahedron_mesh
    TetrahedronMesh.cpp
)

