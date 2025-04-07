include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver2D)
include_directories(${CMAKE_SOURCE_DIR}/include/LinearDGSolver3D)


# 设置目标目录
set(PY_MODULE_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/pymodule/testPyModule)

# 创建目标目录（构建时）
file(MAKE_DIRECTORY ${PY_MODULE_OUTPUT_DIR})

# PY_DG_2D 模块
pybind11_add_module(PY_DG_2D DG_2D_binding.cpp)
target_link_libraries(PY_DG_2D PRIVATE LinearDGSolver2DCycleBoundary LinearDGSolver2DCavityCylinderFlow triangle_mesh GeometryUtils)

add_custom_command(TARGET PY_DG_2D POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PY_DG_2D> ${PY_MODULE_OUTPUT_DIR}
    COMMENT "Copying PY_DG_2D to ${PY_MODULE_OUTPUT_DIR}"
)

# PY_DG_3D 模块
pybind11_add_module(PY_DG_3D DG_3D_binding.cpp)
target_link_libraries(PY_DG_3D PRIVATE LinearDGSolver3DCycleBoundary tetrahedron_mesh GeometryUtils)

add_custom_command(TARGET PY_DG_3D POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PY_DG_3D> ${PY_MODULE_OUTPUT_DIR}
    COMMENT "Copying PY_DG_3D to ${PY_MODULE_OUTPUT_DIR}"
)