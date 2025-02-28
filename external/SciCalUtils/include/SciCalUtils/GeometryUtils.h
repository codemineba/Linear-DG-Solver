#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H

#include <iostream>
#include <cmath>
#include <algorithm>

// 用于判断两个 double 是否近似相等
bool is_almost_equal(double a, double b, double scale = 1e-9);

// 冒泡排序函数：对 ys 数组排序，并同步更新 indices 数组
void bubble_sort_with_indices(double *arr, int *indices, int size);

// 计算两个向量的点积
double dot(double* a, double* b, unsigned long size);

// 计算一个向量自身的L2范数
double l2_norm(double*a, unsigned long size);

// 计算两个向量的欧式距离
double euclidean_distance(double* a, double* b, unsigned long size);

// 计算两个向量的叉积
void cross(const double u[3], const double v[3], double result[3]);

// 计算二维三角形面积 (带方向)
double directed_2D_triangle_area(double* v0, double* v1, double* v2);

// 计算三维三角形面积 (带方向)
double directed_3D_triangle_area(double* v0, double* v1, double* v2);

// 计算三角形关于顶点坐标的方向面积梯度
void directed_2D_triangle_area_gradient(double* v0, double* v1, double* v2, int idx, double* area_v);

// 计算四面体体积 (带方向)
double directed_3D_tetrahedron_volume(double* v0, double* v1, double* v2, double* v3);

// 计算四面体关于顶点坐标的方向体积梯度
void directed_3D_tetrahedron_volume_gradient(double* v0, double* v1, double* v2, double* v3, int idx, double* volume_v);

// 判断两条竖直(水平)边是否沿水平(竖直)方向平移
bool is_edge_shift_from(double edge1[2][2], double edge2[2][2]);

// 判断两三角形是否沿x, y, z方向平移
bool is_tri_shift_from(double tri1[3][3], double tri2[3][3]);

#endif // GEOMETRY_UTILS_H