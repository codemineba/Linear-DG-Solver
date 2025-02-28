include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver2D)
include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver3D)

# 将 CMAKE_SOURCE_DIR 传递到代码中作为宏
add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

# 测试可执行文件
add_executable(2D_LDG_CB 2D_LDG_CB.cpp)
target_link_libraries(2D_LDG_CB LinearDGSolver2DCycleBoundary triangle_mesh GeometryUtils)

add_executable(2D_LDG_CCF 2D_LDG_CCF.cpp)
target_link_libraries(2D_LDG_CCF LinearDGSolver2DCavityCylinderFlow triangle_mesh GeometryUtils)

add_executable(3D_LDG_CB 3D_LDG_CB.cpp)
target_link_libraries(3D_LDG_CB LinearDGSolver3DCycleBoundary tetrahedron_mesh GeometryUtils)

