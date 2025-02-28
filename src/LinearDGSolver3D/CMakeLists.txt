include_directories(${CMAKE_SOURCE_DIR}/include/mesh_structure)
include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver3D)
include_directories(${CMAKE_SOURCE_DIR}/include/Utils)

# 3D_LDG_CB
add_library(LinearDGSolver3DCycleBoundary
    LinearDGSolver3D.cpp
    LinearDGSolver3DCycleBoundary.cpp
)

