include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver2D)

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

# 🔗 链接外部库
target_link_libraries(LinearDGSolver2DCycleBoundary PRIVATE GeometryUtils triangle_mesh)
target_link_libraries(LinearDGSolver2DCavityCylinderFlow PRIVATE GeometryUtils triangle_mesh)