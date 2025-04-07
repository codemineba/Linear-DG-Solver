#include "LinearDGSolver2DCavityObstacleFlow.h"


// 计算边界类型
void LinearDGSolver_2D_CavityObstacleFlow::computeBoundaryType() {
    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    unsigned long *boundary=trimesh_->boundary();
    unsigned long ** edgeidx=trimesh_->edge_info();

    // 判断网格是否为fang'qiang

    boundary_type_ = new unsigned long[nEdge];  
    for(unsigned long i=0; i < nEdge; i++){   // 先初始化所有的边界信息为边数+1 表示内部边
        boundary_type_[i] = nEdge+1;
    }
    
    // 判断边界类型
    for(unsigned long i=0; i < nBoundary; i++){
        unsigned long idx_edge = boundary[i];

        // 该边的两个点
        unsigned long p1 = edgeidx[0][idx_edge];
        unsigned long p2 = edgeidx[1][idx_edge];
        double vx1 = x[p1];
        double vy1 = y[p1];
        double vx2 = x[p2];
        double vy2 = y[p2];
    
        // 判断该边界是与x轴平行还是与y轴平行, 并查找对应的周期边界
        double scale = 1.0e-10;   // 浮点数不能直接用与比较大小, 定义一个小的范围来进行比较
        if(fabs(vx1-vx2) < scale && fabs(vx1-left_boundary_value_) < scale){  // 左边界 (入流边界)
            boundary_type_[idx_edge] = 0;
        }else if(((fabs(vx1-vx2) < scale) && (fabs(vx1-right_boundary_value_)< scale)) || 
         ((fabs(vy1-vy2) < scale) && (fabs(vy1-top_boundary_value_)< scale)) ||
         ((fabs(vy1-vy2) < scale) && (fabs(vy1-bottom_boundary_value_)< scale))){  // 右边界, 上下边界 (流出边界)
            boundary_type_[idx_edge] = 1;  
        }else{  // 障碍物边界 (反射边界)
            boundary_type_[idx_edge] = 2;
        }
    }
}


// 计算单元边界数值流通量
void LinearDGSolver_2D_CavityObstacleFlow::computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[2][4]) {   
    unsigned long** tri_edge_conn = trimesh_->tri_edge_connection();

    // 边编号
    unsigned long idx_edge = tri_edge_conn[ie][ik];

    // 外法向量
    double normal[2];
    normal[0] = outerNormal_[0][ik*3+ie];
    normal[1] = outerNormal_[1][ik*3+ie];

    // 计算u(x, y) 在高斯点上的值
    double u1[nBarx_][4], u2[nBarx_][4];
    // 出流
    computeNumericalSolutionOnGaussPoints(ik, ie, u, u1);

    // 入流
    // 根据不同的边界类型做出不同的处理
    if(boundary_type_[idx_edge] == 0){  // 入流边界
        // 入流为固定的初始条件
        computeNumericalSolutionOnGaussPoints(ik, ie, u0_, u2);  // u0_即使用初始的u值计算
    }else if(boundary_type_[idx_edge] == 1){  // 出流边界
        // 入流等于出流
        std::memcpy(u2, u1, sizeof(u1));  // u2=u1  
    }else if(boundary_type_[idx_edge] == 2){  // 反射边界
        // 密度与能量相等 速度法向相反 切向相同
        for(int i=0; i<nBarx_; i++){
            u2[i][0] = u1[i][0];
            u2[i][3] = u1[i][3];
            u2[i][1] = -normal[0]*(u1[i][1]*normal[0]+u1[i][2]*normal[1]) - normal[1]*(-u1[i][1]*normal[1]+u1[i][2]*normal[0]);
            u2[i][2] = -normal[1]*(u1[i][1]*normal[0]+u1[i][2]*normal[1]) + normal[0]*(-u1[i][1]*normal[1]+u1[i][2]*normal[0]);  
        }
    }   
    
    // 计算数值流通量
    for (int m=0; m<nBarx_; m++){  // 对于两个高斯点
        double ff[nVars_];
        numerical_flux(u1[m], u2[m], normal, ff);
        for(int n=0; n<nVars_; n++){ // 对于4个物理量
            f[m][n] = ff[n]; 
        }
    }
}

// 计算Rho的L2范数
double LinearDGSolver_2D_CavityObstacleFlow::computeL2NormOfRho(){

    double rho[nVertex];

    // 密度数值解
    int idx = 0;
    for(unsigned long i=0; i<nVertex; i++){ rho[i] = vertexFieldValues[i+idx*nVertex]; }

    double L2Result = l2_norm(rho, nVertex)/nVertex;

    return L2Result;
}

// 输出日志信息 (rho误差的L2范数)
std::string LinearDGSolver_2D_CavityObstacleFlow::LogMessage() {
    std::ostringstream oss;
    // 设置输出格式
    oss << "dt is " <<timeInterval_ << ", rho's L2Normal is " << computeL2NormOfRho();
    return oss.str();
}