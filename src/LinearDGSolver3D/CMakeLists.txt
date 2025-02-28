include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver3D)


# 3D_LDG_CB
add_library(LinearDGSolver3DCycleBoundary
    LinearDGSolver3D.cpp
    LinearDGSolver3DCycleBoundary.cpp
)

# 🔗 链接外部库
target_link_libraries(LinearDGSolver3DCycleBoundary PRIVATE GeometryUtils triangle_mesh)