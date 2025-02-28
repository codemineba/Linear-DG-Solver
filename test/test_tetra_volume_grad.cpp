#include <iostream>

// 计算三角形面积 (带方向)
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
    return area;  // 返回带方向的面积
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

int main() {
    // 四面体顶点坐标
    double v0[] = {0.0, 0.0, 0.0};
    double v1[] = {1.0, 0.0, 0.0};
    double v2[] = {0.0, 1.0, 0.0};
    double v3[] = {0.0, 0.0, 1.0};

//    double v0[] = {0.0, 0.0, 0.0};
//    double v1[] = {100.0, 0.0, 0.0};
//    double v2[] = {0.0, 100.0, 0.0};
//    double v3[] = {0.0, 100.0, 1.0};
    
//    double v0[] = {10.0, 10.0, 10.0};
//    double v1[] = {20.0, 10.0, 10.0};
//    double v2[] = {10.0, 20.0, 10.0};
//    double v3[] = {10.0, 10.0, 30.0};

    double volume_v[3]; // 体积梯度数组

    // 测试计算四面体体积及其梯度
    double volume = directed_3D_tetrahedron_volume(v0, v1, v2, v3);
    std::cout << "Tetrahedron Volume: " << volume << std::endl;

    // 计算关于每个顶点的体积梯度
    for (int i = 0; i < 4; i++) {
        directed_3D_tetrahedron_volume_gradient(v0, v1, v2, v3, i, volume_v);
        std::cout << "Gradient w.r.t. vertex " << i << " : (" 
                  << volume_v[0] << ", " 
                  << volume_v[1] << ", " 
                  << volume_v[2] << ")" << std::endl;
    }

    return 0;
}
