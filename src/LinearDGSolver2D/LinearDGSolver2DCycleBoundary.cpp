#include "LinearDGSolver2DCycleBoundary.h"


// 计算周期边界信息
void LinearDGSolver_2D_CycleBoundary::computePeriodicBoundaryInfo(){
    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    unsigned long *boundary = trimesh_->boundary();
    unsigned long ** edgeidx=trimesh_->edge_info();
    unsigned long** edge_order_in_tri = trimesh_->edge_order_in_tri();

    // 周期边界
    periodic_boundary_ = new unsigned long[nEdge];  
    for(unsigned long i=0; i < nEdge; i++){   // 先初始化所有的边界信息为边数+1 表示内部边
        periodic_boundary_[i] = nEdge+1;
    }
    
    // 寻找边界的对应边的索引 （周期边界） 
    // 记录周期边界信息(要求网格点在边界处一致分布)
    for(unsigned long i=0; i < nBoundary; i++){
        // 该边的两个点
        unsigned long ie_p1 = edgeidx[0][boundary[i]];
        unsigned long ie_p2 = edgeidx[1][boundary[i]];
        double edge1[2][2]={{x[ie_p1], y[ie_p1]}, {x[ie_p2], y[ie_p2]}};

        for(unsigned long j=0; j < nBoundary; j++){
            // 该边的两个点
            unsigned long je_p1 = edgeidx[0][boundary[j]];
            unsigned long je_p2 = edgeidx[1][boundary[j]];
            double edge2[2][2]={{x[je_p1], y[je_p1]}, {x[je_p2], y[je_p2]}};

            if(is_edge_shift_from(edge1, edge2)){  // 符合平移
                periodic_boundary_[boundary[i]]= boundary[j]; 
                break;
            }
        }
    }
    // 调整周期边界之间高斯点的顺序 (确保计算边界上的数值通量时,高斯点位置能直接对应)
    for(unsigned long i=0; i<nBoundary; i++){
        // 边界上的两个高斯点
        double v1[2] = {barx_[0][nBarx_*boundary[i]], barx_[1][nBarx_*boundary[i]]};
        double v2[2] = {barx_[0][nBarx_*boundary[i]+1], barx_[1][nBarx_*boundary[i]+1]};
        if((is_almost_equal(v1[0], v2[0]) && (v1[1] > v2[1])) || // 关于y轴平行 规定高斯点从下往上排列(y轴正方向) 即y0 < y1
        (is_almost_equal(v1[1], v2[1]) && (v1[0] > v2[0])))      // 关于x轴平行 规定高斯点从左往右排列(x轴正方向) 即x0 < x1
        {  
            // 交换 v1 和 v2 的坐标 以及基函数的值
            unsigned long ik = edgeidx[2][boundary[i]]; 
            unsigned long order = edge_order_in_tri[0][boundary[i]];
            std::swap(barx_[0][nBarx_ * boundary[i]], barx_[0][nBarx_ * boundary[i] + 1]);
            std::swap(barx_[1][nBarx_ * boundary[i]], barx_[1][nBarx_ * boundary[i] + 1]);
            std::swap(phi_barx_[0][nBarx_*(nSide_*ik+order)], phi_barx_[0][nBarx_*(nSide_*ik+order)+1]);  
            std::swap(phi_barx_[1][nBarx_*(nSide_*ik+order)], phi_barx_[1][nBarx_*(nSide_*ik+order)+1]);  
            std::swap(phi_barx_[2][nBarx_*(nSide_*ik+order)], phi_barx_[2][nBarx_*(nSide_*ik+order)+1]);  
            
        }
    }
}


// 计算初值条件
void LinearDGSolver_2D_CycleBoundary::computeInitialCondition(unsigned long n, double *x, double *y, double **rtval) {
    for (unsigned long i=0; i<n; i++){
        // 初始条件
        rtval[0][i] = 1+0.2*sin(PI*(x[i]+y[i]));
        rtval[1][i] = 0.7*rtval[0][i];
        rtval[2][i] = 0.3*rtval[0][i];
        rtval[3][i] = 1.0/(GAMMA-1)+1.0/2.0*rtval[0][i]*(0.7*0.7+0.3*0.3);
    }
}


// 计算单元边界数值流通量
void LinearDGSolver_2D_CycleBoundary::computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[2][4]) {   
    unsigned long** edgeInfo = trimesh_->edge_info();
    unsigned long** tri_edge_conn = trimesh_->tri_edge_connection();
    unsigned long** edge_order_in_tri = trimesh_->edge_order_in_tri();

    // 边编号
    unsigned long idx_edge = tri_edge_conn[ie][ik];

    // 外法向量
    double normal[2];
    normal[0] = outerNormal_[0][ik*nSide_+ie];
    normal[1] = outerNormal_[1][ik*nSide_+ie];

    // 周期边界
    unsigned long idx_edge_ = periodic_boundary_[idx_edge];

    if (idx_edge_ == nEdge+1)  // 如果是内部边则抛出异常
        std::cerr << "Error: could not find the corresponding boundary" << std::endl; 

    // idx_edge_所在的单元 (由于idx_edge_也是边界 所以只存在一个与其相连的单元)
    unsigned ik_ = edgeInfo[2][idx_edge_];
    // idx_edge_在ik_中的顺序ie_
    unsigned long ie_ = edge_order_in_tri[0][idx_edge_];


    // 计算u(x, y) 在高斯点上的值
    double u1[nBarx_][4], u2[nBarx_][4];
    computeNumericalSolutionOnGaussPoints(ik, ie, u, u1);
    computeNumericalSolutionOnGaussPoints(ik_, ie_, u, u2);
    
    // 计算数值流通量
    for (int m=0; m<nBarx_; m++){  // 对于两个高斯点
        double ff[nVars_];
        numerical_flux(u1[m], u2[m], normal, ff);
        for(int n=0; n<nVars_; n++){ // 对于4个物理量
            f[m][n] = ff[n]; 
        }
    }
}


// 误差估计
double LinearDGSolver_2D_CycleBoundary::computeL2ErrorOfRho(){
    
    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    double rho[nVertex];

    // 密度数值解
    int idx = 0;
    for(unsigned long i=0; i<nVertex; i++){ rho[i] = vertexFieldValues[i+idx*nVertex]; }

    // 密度准确解
    double exactSolution[nVertex];
    for (unsigned long i =0; i<nVertex; i++){
        exactSolution[i] = 1.0 + 0.2 * sin(PI*(x[i] + y[i] - timePoint_)); 
    }

    // 计算误差的平均l2范数
    double L2Err = euclidean_distance(exactSolution, rho, nVertex)/nVertex;

    return L2Err;
}

// 输出日志信息 (rho误差的L2范数)
std::string LinearDGSolver_2D_CycleBoundary::LogMessage() {
    std::ostringstream oss;
    // 设置输出格式
    oss << "dt is " <<timeInterval_ << ", rho's L2Error is " << computeL2ErrorOfRho();
    return oss.str();
}
