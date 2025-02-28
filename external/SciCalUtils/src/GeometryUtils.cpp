#include "SciCalUtils/GeometryUtils.h"

// 定义一个用于判断两个 double 是否近似相等的函数
bool is_almost_equal(double a, double b, double scale) {
    return abs(a - b) < scale;
}


// 冒泡排序函数：对 ys 数组排序，并同步更新 indices 数组
void bubble_sort_with_indices(double *arr, int *indices, int size) {
    for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size - 1 - i; ++j) {
            if (arr[j] > arr[j + 1]) {
                std::swap(arr[j], arr[j + 1]);
                std::swap(indices[j], indices[j + 1]);
            }
        }
    }
}


// 计算两个向量的点积
double dot(double* a, double* b, unsigned long size) {
    double result = 0;
    for (unsigned long i = 0; i < size; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// 计算一个向量自身的L2范数
double l2_norm(double*a, unsigned long size){
    return sqrt(dot(a, a, size));
}

// 计算两个向量的欧式距离
double euclidean_distance(double* a, double* b, unsigned long size){
    double c[size];
    for(unsigned long i=0; i<size; i++){
        c[i] = a[i]-b[i];
    }
    return l2_norm(c, size);
}

// 计算两个向量的叉积 (默认参数三维向量 返回一个三维向量)
void cross(const double u[3], const double v[3], double result[3]){
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
}

// 计算二维三角形面积 (带方向)
double directed_2D_triangle_area(double* v0, double* v1, double* v2) {
    double x0 = v0[0], y0 = v0[1];
    double x1 = v1[0], y1 = v1[1];
    double x2 = v2[0], y2 = v2[1];

    // 方向性由符号体现 逆时针为正
    double area = 0.5 * (
        x0 * (y1 - y2) +
        x1 * (y2 - y0) +
        x2 * (y0 - y1)
    );
    return area; // 返回带方向的面积
}

// 计算三维三角形面积 (带方向)
double directed_3D_triangle_area(double* v0, double* v1, double* v2) {
    // 构造向量
    double AB[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
    double AC[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

    // 叉积的模长的一半即为空间三角形的面积
    double normal[3]; 
    cross(AB, AC, normal);
    double area = l2_norm(normal, 3)/2;

    return area; // 返回带方向的面积
}

// 计算关于顶点坐标的方向面积梯度
void directed_2D_triangle_area_gradient(double* v0, double* v1, double* v2, int idx, double* area_v){
    double x0 = v0[0], y0 = v0[1];
    double x1 = v1[0], y1 = v1[1];
    double x2 = v2[0], y2 = v2[1];

    if (idx==0){   // 关于v_idx的面积梯度
        area_v[0]=0.5*(y1-y2), area_v[1]=0.5*(x2-x1);  // dA/dv0
    }else if(idx==1){
        area_v[0]=0.5*(y2-y0), area_v[1]=0.5*(x0-x2);  // dA/dv1
    }else if(idx==2){
        area_v[0]=0.5*(y0-y1), area_v[1]=0.5*(x1-x0);  // dA/dv2
    }else{
        std::cerr<<"gradient of area error! please enter the correct idx 0,1,2"<<std::endl;
        throw -1;
    }
}

// 计算四面体体积 (带方向)
double directed_3D_tetrahedron_volume(double* v0, double* v1, double* v2, double* v3) {
    double v0_2d[2]={v0[0], v0[1]};
    double v1_2d[2]={v1[0], v1[1]};
    double v2_2d[2]={v2[0], v2[1]};
    double v3_2d[2]={v3[0], v3[1]};
    
    double z0=v0[2], z1=v1[2], z2=v2[2], z3=v3[2]; 

    double volume = (1.0 / 3.0) * (
        -1 * z0 * directed_2D_triangle_area(v1_2d, v2_2d, v3_2d)
        +1 * z1 * directed_2D_triangle_area(v0_2d, v2_2d, v3_2d)
        -1 * z2 * directed_2D_triangle_area(v0_2d, v1_2d, v3_2d)
        +1 * z3 * directed_2D_triangle_area(v0_2d, v1_2d, v2_2d)
    );
    return volume;
}


// 计算关于顶点坐标的方向体积梯度
void directed_3D_tetrahedron_volume_gradient(double* v0, double* v1, double* v2, double* v3, int idx, double* volume_v){
    double v0_2d[2]={v0[0], v0[1]};
    double v1_2d[2]={v1[0], v1[1]};
    double v2_2d[2]={v2[0], v2[1]};
    double v3_2d[2]={v3[0], v3[1]};
    
    double z0=v0[2], z1=v1[2], z2=v2[2], z3=v3[2]; 

    if (idx==0){   // 关于v_idx的体积梯度
    	double s1_v0[2], s2_v0[2], s3_v0[2];  
    	directed_2D_triangle_area_gradient(v0_2d, v2_2d, v3_2d, 0, s1_v0);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v3_2d, 0, s2_v0);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v2_2d, 0, s3_v0);
    	// dV/dv0
        volume_v[0]=(1.0/3.0) * (z1*s1_v0[0]-z2*s2_v0[0]+z3*s3_v0[0]);
		volume_v[1]=(1.0/3.0) * (z1*s1_v0[1]-z2*s2_v0[1]+z3*s3_v0[1]);  
		volume_v[2]=(1.0/3.0) * (-directed_2D_triangle_area(v1_2d, v2_2d, v3_2d));  
    }else if(idx==1){
        double s0_v1[2], s2_v1[2], s3_v1[2];  
    	directed_2D_triangle_area_gradient(v1_2d, v2_2d, v3_2d, 0, s0_v1);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v3_2d, 1, s2_v1);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v2_2d, 1, s3_v1);
    	// dV/dv1
        volume_v[0]=(1.0/3.0) * (-z0*s0_v1[0]-z2*s2_v1[0]+z3*s3_v1[0]);
		volume_v[1]=(1.0/3.0) * (-z0*s0_v1[1]-z2*s2_v1[1]+z3*s3_v1[1]);  
		volume_v[2]=(1.0/3.0) * (directed_2D_triangle_area(v0_2d, v2_2d, v3_2d));  
    }else if(idx==2){
        double s0_v2[2], s1_v2[2], s3_v2[2];  
    	directed_2D_triangle_area_gradient(v1_2d, v2_2d, v3_2d, 1, s0_v2);
    	directed_2D_triangle_area_gradient(v0_2d, v2_2d, v3_2d, 1, s1_v2);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v2_2d, 2, s3_v2);
    	// dV/dv2
        volume_v[0]=(1.0/3.0) * (-z0*s0_v2[0]+z1*s1_v2[0]+z3*s3_v2[0]);
		volume_v[1]=(1.0/3.0) * (-z0*s0_v2[1]+z1*s1_v2[1]+z3*s3_v2[1]);  
		volume_v[2]=(1.0/3.0) * (-directed_2D_triangle_area(v0_2d, v1_2d, v3_2d));  
	}else if(idx==3){
        double s0_v3[2], s1_v3[2], s2_v3[2];  
    	directed_2D_triangle_area_gradient(v1_2d, v2_2d, v3_2d, 2, s0_v3);
    	directed_2D_triangle_area_gradient(v0_2d, v2_2d, v3_2d, 2, s1_v3);
    	directed_2D_triangle_area_gradient(v0_2d, v1_2d, v3_2d, 2, s2_v3);
    	// dV/dv3
        volume_v[0]=(1.0/3.0) * (-z0*s0_v3[0]+z1*s1_v3[0]-z2*s2_v3[0]);
		volume_v[1]=(1.0/3.0) * (-z0*s0_v3[1]+z1*s1_v3[1]-z2*s2_v3[1]);  
		volume_v[2]=(1.0/3.0) * (directed_2D_triangle_area(v0_2d, v1_2d, v2_2d));  
    }else{
        std::cerr<<"gradient of area error! please enter the correct idx 0,1,2,3"<<std::endl;
        throw -1;
    }
}

// 判断 edge1 是否可以通过平移得到 edge2
bool is_edge_shift_from(double edge1[2][2], double edge2[2][2]) {
    double edge1_vx1 = edge1[0][0], edge1_vy1 = edge1[0][1], edge1_vx2 = edge1[1][0], edge1_vy2 = edge1[1][1];
    double edge2_vx1 = edge2[0][0], edge2_vy1 = edge2[0][1], edge2_vx2 = edge2[1][0], edge2_vy2 = edge2[1][1]; 
    // 检查 edge1 和 edge2 是否为水平边
    bool is_edge1_horizontal = is_almost_equal(edge1_vy1, edge1_vy2);
    bool is_edge2_horizontal = is_almost_equal(edge2_vy1, edge2_vy2);

    // 检查 edge1 和 edge2 是否为竖直边
    bool is_edge1_vertical = is_almost_equal(edge1_vx1, edge1_vx2);
    bool is_edge2_vertical = is_almost_equal(edge2_vx1, edge2_vx2);

    // 如果 edge1 和 edge2 同为水平边
    if (is_edge1_horizontal && is_edge2_horizontal) {
        // 检查两个点的 x 坐标是否完全相等且不重合
        return  (is_almost_equal(edge2_vx1, edge1_vx1)||is_almost_equal(edge2_vx1, edge1_vx2)) &&
         (is_almost_equal(edge2_vx2, edge1_vx1)||is_almost_equal(edge2_vx2, edge1_vx2)) &&
          !(is_almost_equal(edge1_vy1, edge2_vy1));
    }

    // 如果 edge1 和 edge2 同为竖直边
    if (is_edge1_vertical && is_edge2_vertical) {
        // 检查两个点的 y 坐标是否完全相等且不重合
        return  (is_almost_equal(edge2_vy1, edge1_vy1)||is_almost_equal(edge2_vy1, edge1_vy2)) &&
         (is_almost_equal(edge2_vy2, edge1_vy1)||is_almost_equal(edge2_vy2, edge1_vy2)) &&
          !(is_almost_equal(edge1_vx1, edge2_vx1));
    }

    // 不满足上述条件
    return false;
}


// 判断 tri1 是否可以通过平移得到 tri2
bool is_tri_shift_from(double tri1[3][3], double tri2[3][3]) {
    double tri1_vx1 = tri1[0][0], tri1_vy1 = tri1[0][1], tri1_vz1 = tri1[0][2];
    double tri1_vx2 = tri1[1][0], tri1_vy2 = tri1[1][1], tri1_vz2 = tri1[1][2];
    double tri1_vx3 = tri1[2][0], tri1_vy3 = tri1[2][1], tri1_vz3 = tri1[2][2];

    double tri2_vx1 = tri2[0][0], tri2_vy1 = tri2[0][1], tri2_vz1 = tri2[0][2];
    double tri2_vx2 = tri2[1][0], tri2_vy2 = tri2[1][1], tri2_vz2 = tri2[1][2];
    double tri2_vx3 = tri2[2][0], tri2_vy3 = tri2[2][1], tri2_vz3 = tri2[2][2];

    // 检查 tri1 和 tri2 是否平行于 xoy 平面
    bool is_tri1_parallel_to_xoy = is_almost_equal(tri1_vz1, tri1_vz2) && is_almost_equal(tri1_vz1, tri1_vz3);
    bool is_tri2_parallel_to_xoy = is_almost_equal(tri2_vz1, tri2_vz2) && is_almost_equal(tri2_vz1, tri2_vz3);

    // 检查 tri1 和 tri2 是否平行于 yoz 平面
    bool is_tri1_parallel_to_yoz = is_almost_equal(tri1_vx1, tri1_vx2) && is_almost_equal(tri1_vx1, tri1_vx3);
    bool is_tri2_parallel_to_yoz = is_almost_equal(tri2_vx1, tri2_vx2) && is_almost_equal(tri2_vx1, tri2_vx3);

    // 检查 tri1 和 tri2 是否平行于 zox 平面
    bool is_tri1_parallel_to_zox = is_almost_equal(tri1_vy1, tri1_vy2) && is_almost_equal(tri1_vy1, tri1_vy3);
    bool is_tri2_parallel_to_zox = is_almost_equal(tri2_vy1, tri2_vy2) && is_almost_equal(tri2_vy1, tri2_vy3);


    // 如果 tri1 和 tri2 同为平行于 xoy 平面 
    if (is_tri1_parallel_to_xoy && is_tri2_parallel_to_xoy) {

        double tri1_vx[3] = {tri1_vx1, tri1_vx2, tri1_vx3};
        double tri1_vy[3] = {tri1_vy1, tri1_vy2, tri1_vy3};
        double tri2_vx[3] = {tri2_vx1, tri2_vx2, tri2_vx3};
        double tri2_vy[3] = {tri2_vy1, tri2_vy2, tri2_vy3};
        std::sort(tri1_vx, tri1_vx + 3);
        std::sort(tri1_vy, tri1_vy + 3);
        std::sort(tri2_vx, tri2_vx + 3);
        std::sort(tri2_vy, tri2_vy + 3);

        // 检查两个三角形的 x 和 y 坐标是否完全相等且不重合
        return (is_almost_equal(tri2_vx[0], tri1_vx[0]) && is_almost_equal(tri2_vx[1], tri1_vx[1]) && is_almost_equal(tri2_vx[2], tri1_vx[2]) &&
            is_almost_equal(tri2_vy[0], tri1_vy[0]) && is_almost_equal(tri2_vy[1], tri1_vy[1]) && is_almost_equal(tri2_vy[2], tri1_vy[2]) &&
            !(is_almost_equal(tri1_vz1, tri2_vz1)));
    }

    // 如果 tri1 和 tri2 同为平行于 yoz 平面 
    if (is_tri1_parallel_to_yoz && is_tri2_parallel_to_yoz) {
        double tri1_vy[3] = {tri1_vy1, tri1_vy2, tri1_vy3};
        double tri1_vz[3] = {tri1_vz1, tri1_vz2, tri1_vz3};
        double tri2_vy[3] = {tri2_vy1, tri2_vy2, tri2_vy3};
        double tri2_vz[3] = {tri2_vz1, tri2_vz2, tri2_vz3};
        
        std::sort(tri1_vy, tri1_vy + 3);
        std::sort(tri1_vz, tri1_vz + 3);
        std::sort(tri2_vy, tri2_vy + 3);
        std::sort(tri2_vz, tri2_vz + 3);

        // 检查两个三角形的 y 和 z 坐标是否完全相等且不重合
        return (is_almost_equal(tri2_vy[0], tri1_vy[0]) && is_almost_equal(tri2_vy[1], tri1_vy[1]) && is_almost_equal(tri2_vy[2], tri1_vy[2]) &&
                is_almost_equal(tri2_vz[0], tri1_vz[0]) && is_almost_equal(tri2_vz[1], tri1_vz[1]) && is_almost_equal(tri2_vz[2], tri1_vz[2]) &&
                !(is_almost_equal(tri1_vx1, tri2_vx1)));
    }

    // 如果 tri1 和 tri2 同为平行于 zox 平面 
    if (is_tri1_parallel_to_zox && is_tri2_parallel_to_zox) {
        double tri1_vz[3] = {tri1_vz1, tri1_vz2, tri1_vz3};
        double tri1_vx[3] = {tri1_vx1, tri1_vx2, tri1_vx3};
        double tri2_vz[3] = {tri2_vz1, tri2_vz2, tri2_vz3};
        double tri2_vx[3] = {tri2_vx1, tri2_vx2, tri2_vx3};

        std::sort(tri1_vz, tri1_vz + 3);
        std::sort(tri1_vx, tri1_vx + 3);
        std::sort(tri2_vz, tri2_vz + 3);
        std::sort(tri2_vx, tri2_vx + 3);

        // 检查两个三角形的 z 和 x 坐标是否完全相等且不重合
        return (is_almost_equal(tri2_vz[0], tri1_vz[0]) && is_almost_equal(tri2_vz[1], tri1_vz[1]) && is_almost_equal(tri2_vz[2], tri1_vz[2]) &&
                is_almost_equal(tri2_vx[0], tri1_vx[0]) && is_almost_equal(tri2_vx[1], tri1_vx[1]) && is_almost_equal(tri2_vx[2], tri1_vx[2]) &&
                !(is_almost_equal(tri1_vy1, tri2_vy1)));
    }


    // 不满足上述条件
    return false;
}