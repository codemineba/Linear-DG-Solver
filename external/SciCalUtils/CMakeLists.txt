cmake_minimum_required(VERSION 3.10)
project(SciCalUtils)

# 创建 GeometryUtils 库
add_library(GeometryUtils src/GeometryUtils.cpp)
target_include_directories(GeometryUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# # 创建 CalUtils 库
# add_library(CalUtils src/CalUtils.cpp)
# target_include_directories(CalUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)