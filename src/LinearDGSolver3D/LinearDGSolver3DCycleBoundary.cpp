#include "LinearDGSolver3DCycleBoundary.h"


// 计算周期边界信息
void LinearDGSolver_3D_CycleBoundary::computePeriodicBoundaryInfo(){
    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();
    unsigned long ** faceidx=tetmesh_->face_info();
    unsigned long** face_order_in_tet = tetmesh_->face_order_in_tet();

    // 周期边界
    periodic_boundary_ = new unsigned long[nFace];  
    for(unsigned long i=0; i < nFace; i++){   // 先初始化所有的边界信息为面数+1 表示内部面
        periodic_boundary_[i] = nFace+1;
    }
    
    // 寻找边界的对应面的索引 （周期边界） 
    // 记录周期边界信息(要求网格点在边界处一致分布)
    for(unsigned long i=0; i < nBoundry; i++){
        // 该三角形的三个点
        unsigned long itri_p1 = faceidx[0][boundary_[i]];
        unsigned long itri_p2 = faceidx[1][boundary_[i]];
        unsigned long itri_p3 = faceidx[2][boundary_[i]];
        double tri1[3][3]={{x[itri_p1], y[itri_p1], z[itri_p1]}, {x[itri_p2], y[itri_p2], z[itri_p2]}, {x[itri_p3], y[itri_p3], z[itri_p3]}};

        for(unsigned long j=0; j < nBoundry; j++){
            // 该边的两个点
            unsigned long jtri_p1 = faceidx[0][boundary_[j]];
            unsigned long jtri_p2 = faceidx[1][boundary_[j]];
            unsigned long jtri_p3 = faceidx[2][boundary_[j]];
            double tri2[3][3]={{x[jtri_p1], y[jtri_p1], z[jtri_p1]}, {x[jtri_p2], y[jtri_p2], z[jtri_p2]}, {x[jtri_p3], y[jtri_p3], z[jtri_p3]}};

            if(is_tri_shift_from(tri1, tri2)){  // 符合平移
                periodic_boundary_[boundary_[i]]= boundary_[j]; 
                break;
            }
        }
    }
    // 调整周期边界之间边界积分点的顺序 (确保计算边界上的数值通量时,边界求积点位置能直接对应)
    for(unsigned long i=0; i<nBoundry; i++){
        // 边界上的三个积分点
        double v1[dim_] = {barx_[0][nBarx_*boundary_[i]], barx_[1][nBarx_*boundary_[i]], barx_[2][nBarx_*boundary_[i]]};
        double v2[dim_] = {barx_[0][nBarx_*boundary_[i]+1], barx_[1][nBarx_*boundary_[i]+1], barx_[2][nBarx_*boundary_[i]+1]};
        double v3[dim_] = {barx_[0][nBarx_*boundary_[i]+2], barx_[1][nBarx_*boundary_[i]+2], barx_[2][nBarx_*boundary_[i]+2]};


        int idx = -1;  // 初始化索引
        if ((is_almost_equal(v1[0], v2[0]) && !(v1[1] < v2[1] && v2[1] < v3[1]))) {
            idx = 1;  // 关于yoz平面平行 规定积分点沿y轴正方向排列 即y0 < y1 < y2
        } else if ((is_almost_equal(v1[1], v2[1]) && !(v1[2] < v2[2] && v2[2] < v3[2]))) {
            idx = 2;  // 关于zox平面平行 规定积分点沿z轴正方向排列 即z0 < z1 < z2
        } else if ((is_almost_equal(v1[2], v2[2]) && !(v1[0] < v2[0] && v2[0] < v3[0]))) {
            idx = 0;  // 关于xoy平面平行 规定积分点沿x轴正方向排列 即x0 < x1 < x2 
        }

        if (idx != -1) {  
            // 调整 v坐标位置 以及基函数的值 

            int indices[] = {0, 1, 2};
            double arr[] = {v1[idx], v2[idx], v3[idx]};
            // 调用排序函数，对 arr 进行排序，同时更新 indices
            bubble_sort_with_indices(arr, indices, 3);

            // 根据排序后的 indices 数组，重新组织对应的数据
            unsigned long ik = faceidx[3][boundary_[i]];
            unsigned long order = face_order_in_tet[0][boundary_[i]];

            // 临时变量存储排序后的 barx_ 和 phi_barx_ 的值
            double temp_barx[dim_][nBarx_];  
            double temp_phi_barx[nPhi_][nBarx_]; 

            // 读取原数组的值到临时变量，并将排序后的值写回原数组
            for (int j = 0; j < nBarx_; ++j) {  // 遍历每个点
                unsigned long original_index = indices[j];  // 排序前的位置索引
                for (int n = 0; n < dim_; ++n) {  
                    temp_barx[n][j] = barx_[n][nBarx_ * boundary_[i] + original_index];
                }
                for (int n = 0; n < nPhi_; ++n) { 
                    temp_phi_barx[n][j] = phi_barx_[n][nBarx_ * (nSide_*ik+order) + original_index];
                }
            }

            // 交换ik面的积分点的顺序
            for (int j = 0; j < nBarx_; ++j) {  // 遍历每个点
                for (int n = 0; n < dim_; ++n) { 
                    barx_[n][nBarx_ * boundary_[i] + j] = temp_barx[n][j];
                }
                for (int n = 0; n < nPhi_; ++n) {
                    phi_barx_[n][nBarx_ * (nSide_*ik+order) + j] = temp_phi_barx[n][j];
                }
            }
        }
    }
}


// 计算初值条件
void LinearDGSolver_3D_CycleBoundary::computeInitialCondition(unsigned long n, double *x, double *y, double *z, double **rtval) {
    double u0[nVars_][dof_];   // n=dof_
    for (unsigned long i=0; i<n; i++){
        // 初始条件
        u0[0][i] = 1+0.2*sin(PI*(x[i]+y[i]+z[i]));
        u0[1][i] = 0.4*u0[0][i];
        u0[2][i] = 0.3*u0[0][i];
        u0[3][i] = 0.3*u0[0][i];
        u0[4][i] = 1.0/(GAMMA-1)+1.0/2.0*u0[0][i]*(0.4*0.4+0.3*0.3+0.3*0.3);
    }
    
    for (unsigned long ik=0; ik<nElement; ik++){
        for(int j=0; j<nPhi_; j++){
            double phiProduct=0.0;  // 基函数内积
            if(j==0){phiProduct=1.0/6.0;}
            else if(j==1){phiProduct=1.0/160.0;}
            else if(j==2){phiProduct=1.0/180.0;}
            else if(j==3){phiProduct=1.0/240.0;}
            double aj = 6.0*volume_[ik]*phiProduct;
            double scale = volume_[ik]/(nHatx_*aj); 
            for(int i=0; i<nVars_; i++){
                rtval[i][ik*nPhi_+j] = 0.0;
                for(int m=0; m<nHatx_; m++){
                    rtval[i][ik*nPhi_+j]+= scale * u0[i][ik*nHatx_+m] *phi_hatx_[j][ik*nHatx_+m];
                }
            }
        }
    }
}


// 计算单元边界数值流通量
void LinearDGSolver_3D_CycleBoundary::computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[3][5]) {   
    unsigned long** faceInfo = tetmesh_->face_info();
    unsigned long** tet_face_conn = tetmesh_->tet_face_connection();
    unsigned long** face_order_in_tet = tetmesh_->face_order_in_tet();

    // 边编号
    unsigned long idx_edge = tet_face_conn[ie][ik];

    // 外法向量
    double normal[dim_];
    normal[0] = outerNormal_[0][ik*nSide_+ie];
    normal[1] = outerNormal_[1][ik*nSide_+ie];
    normal[2] = outerNormal_[2][ik*nSide_+ie];

    // 周期边界
    unsigned long idx_edge_ = periodic_boundary_[idx_edge];
    
    if (idx_edge_ == nFace+1)  // 如果没有则抛出异常
        std::cerr << "Error: could not find the corresponding boundary" << std::endl; 
    // idx_edge_所在的单元 (由于idx_edge_也是边界 所以只存在一个与其相连的单元)
    unsigned ik_ = faceInfo[3][idx_edge_];
    // idx_edge_在ik_中的顺序ie_
    unsigned long ie_ = face_order_in_tet[0][idx_edge_];

    // std::cout<<std::endl;
    // std::cout<<barx_[0][3*(idx_edge)]<<" "<<barx_[1][3*(idx_edge)]<<" "<<barx_[2][3*(idx_edge)]<<std::endl;
    // std::cout<<barx_[0][3*(idx_edge_)]<<" "<<barx_[1][3*(idx_edge_)]<<" "<<barx_[2][3*(idx_edge_)]<<std::endl;

    // 计算u(x, y, z) 在高斯点上的值
    double u1[nBarx_][5], u2[nBarx_][5];
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
double LinearDGSolver_3D_CycleBoundary::computeL2ErrorOfRho(){
    
    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();

    // 密度准确解
    double* exactSolution = new double[nVertex];
    for (unsigned long i =0; i<nVertex; i++){
        exactSolution[i] = 1.0 + 0.2 * sin(PI*(x[i] + y[i] + z[i] - timePoint_)); 
    }

    // 计算误差的平均l2范数
    double L2Err = euclidean_distance(exactSolution, rho, nVertex)/nVertex;

    return L2Err;
}

// 输出日志信息 (rho误差的L2范数)
std::string LinearDGSolver_3D_CycleBoundary::LogMessage() {
    std::ostringstream oss;
    // 设置输出格式
    oss << "dt is " <<timeInterval_ << ", rho's L2Error is " << computeL2ErrorOfRho();
    return oss.str();
}
