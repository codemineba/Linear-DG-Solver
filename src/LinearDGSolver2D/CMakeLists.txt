include_directories(${CMAKE_SOURCE_DIR}/include/mesh_structure)
include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver2D)
include_directories(${CMAKE_SOURCE_DIR}/include/Utils)

# 2D_LDG_CB
add_library(LinearDGSolver2DCycleBoundary
    LinearDGSolver2D.cpp
    LinearDGSolver2DCycleBoundary.cpp
)

# 2D_LDG_CCF
add_library(LinearDGSolver2DCavityCylinderFlow
    LinearDGSolver2D.cpp
    LinearDGSolver2DCavityObstacleFlow.cpp
    LinearDGSolver2DCavityCylinderFlow.cpp
)
