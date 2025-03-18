#include "LinearDGSolver3DCycleBoundary.h"


// 计算周期边界信息
void LinearDGSolver_3D_CycleBoundary::computePeriodicBoundaryInfo(){
    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();
    unsigned long *boundary=tetmesh_->boundary();
    unsigned long ** faceidx=tetmesh_->face_info();
    unsigned long** face_order_in_tet = tetmesh_->face_order_in_tet();

    // 周期边界
    periodic_boundary_ = new unsigned long[nFace];  
    for(unsigned long i=0; i < nFace; i++){   // 先初始化所有的边界信息为面数+1 表示内部面
        periodic_boundary_[i] = nFace+1;
    }
    
    // 寻找边界的对应面的索引 （周期边界） 
    // 记录周期边界信息(要求网格点在边界处一致分布)
    for(unsigned long i=0; i < nBoundary; i++){
        // 该三角形的三个点
        unsigned long itri_p1 = faceidx[0][boundary[i]];
        unsigned long itri_p2 = faceidx[1][boundary[i]];
        unsigned long itri_p3 = faceidx[2][boundary[i]];
        double tri1[3][3]={{x[itri_p1], y[itri_p1], z[itri_p1]}, {x[itri_p2], y[itri_p2], z[itri_p2]}, {x[itri_p3], y[itri_p3], z[itri_p3]}};

        for(unsigned long j=0; j < nBoundary; j++){
            // 该边的两个点
            unsigned long jtri_p1 = faceidx[0][boundary[j]];
            unsigned long jtri_p2 = faceidx[1][boundary[j]];
            unsigned long jtri_p3 = faceidx[2][boundary[j]];
            double tri2[3][3]={{x[jtri_p1], y[jtri_p1], z[jtri_p1]}, {x[jtri_p2], y[jtri_p2], z[jtri_p2]}, {x[jtri_p3], y[jtri_p3], z[jtri_p3]}};

            if(is_tri_shift_from(tri1, tri2)){  // 符合平移
                periodic_boundary_[boundary[i]]= boundary[j]; 
                break;
            }
        }
    }

    // 调整周期边界之间边界积分点的顺序 (确保计算边界上的数值通量时,边界求积点位置能直接对应)
    double used_bdy[nBoundary];
    int used_times=0;
    for(unsigned long i=0; i<nBoundary; i++){ used_bdy[i]=nFace+1;}
    for(unsigned long bdy_idx=0; bdy_idx<nBoundary; bdy_idx++){

        // 检查是否已经遍历过该边界
        unsigned long bdy = boundary[bdy_idx];
        unsigned long bdy_ = periodic_boundary_[bdy];  // 对应面的边界
        bool has_used=false;
        for(unsigned long i=0; i<nBoundary; i++)
            if(bdy==used_bdy[i]) has_used=true;
        if(has_used) break; 
        used_bdy[used_times*2]=bdy; used_bdy[used_times*2+1]=bdy_;
        used_times++;

        // 边界上的三个积分点
        double v1[3] = {barx_[0][nBarx_*bdy], barx_[1][nBarx_*bdy], barx_[2][nBarx_*bdy]};
        double v2[3] = {barx_[0][nBarx_*bdy+1], barx_[1][nBarx_*bdy+1], barx_[2][nBarx_*bdy+1]};
        double v3[3] = {barx_[0][nBarx_*bdy+2], barx_[1][nBarx_*bdy+2], barx_[2][nBarx_*bdy+2]};
        // 对应边界上的三个积分点
        double v1_[3] = {barx_[0][nBarx_*bdy_], barx_[1][nBarx_*bdy_], barx_[2][nBarx_*bdy_]};
        double v2_[3] = {barx_[0][nBarx_*bdy_+1], barx_[1][nBarx_*bdy_+1], barx_[2][nBarx_*bdy_+1]};
        double v3_[3] = {barx_[0][nBarx_*bdy_+2], barx_[1][nBarx_*bdy_+2], barx_[2][nBarx_*bdy_+2]};

        
        int idx1 = -1, idx2 = -1;   // 初始化索引
        if (are_almost_equal(v1[0], v2[0], v3[0])) {
            idx1=1, idx2=2;  // 关于yoz平面平行 
        } else if (are_almost_equal(v1[1], v2[1], v3[1])) {
            idx1=0, idx2=2;  // 关于zox平面平行 
        } else if (are_almost_equal(v1[2], v2[2], v3[2])) {
            idx1=0, idx2=1;  // 关于xoy平面平行 
        }

        // bdy_的三个积分点按照bdy的顺序排列
        // 共有 3! =6 种可能的排列方式
        int indices[nBarx_];
        double v1_2D[2]={v1[idx1], v1[idx2]}, v2_2D[2]={v2[idx1], v2[idx2]}, v3_2D[2]={v3[idx1], v3[idx2]};
        double v1_2D_[2]={v1_[idx1], v1_[idx2]}, v2_2D_[2]={v2_[idx1], v2_[idx2]}, v3_2D_[2]={v3_[idx1], v3_[idx2]};
        bool adjust=true;
        if (is_almost_equal(v1_2D_, v1_2D, 2) && is_almost_equal(v2_2D_, v2_2D, 2) && is_almost_equal(v3_2D_, v3_2D, 2)) {
            indices[0] = 0; indices[1] = 1; indices[2] = 2; adjust=false;  // 该情况无需调整
        } 
        else if (is_almost_equal(v1_2D_, v1_2D, 2) && is_almost_equal(v2_2D_, v3_2D, 2) && is_almost_equal(v3_2D_, v2_2D, 2)) {
            indices[0] = 0; indices[1] = 2; indices[2] = 1;
        } 
        else if (is_almost_equal(v1_2D_, v2_2D, 2) && is_almost_equal(v2_2D_, v1_2D, 2) && is_almost_equal(v3_2D_, v3_2D, 2)) {
            indices[0] = 1; indices[1] = 0; indices[2] = 2;
        } 
        else if (is_almost_equal(v1_2D_, v2_2D, 2) && is_almost_equal(v2_2D_, v3_2D, 2) && is_almost_equal(v3_2D_, v1_2D, 2)) {
            indices[0] = 1; indices[1] = 2; indices[2] = 0;
        } 
        else if (is_almost_equal(v1_2D_, v3_2D, 2) && is_almost_equal(v2_2D_, v1_2D, 2) && is_almost_equal(v3_2D_, v2_2D, 2)) {
            indices[0] = 2; indices[1] = 0; indices[2] = 1;
        } 
        else if (is_almost_equal(v1_2D_, v3_2D, 2) && is_almost_equal(v2_2D_, v2_2D, 2) && is_almost_equal(v3_2D_, v1_2D, 2)) {
            indices[0] = 2; indices[1] = 1; indices[2] = 0;
        } 
        else {
            std::cerr << "Error: No valid matching found!" << std::endl;
            throw -1;
        }
        
        if (adjust) {  
            std::cout<< "need adjust position of BarX !"<<std::endl;
            // 调整 bdy_ 的坐标位置 以及基函数的值 
            // 根据排序后的 indices 数组，重新组织对应的数据
            unsigned long ik = faceidx[3][bdy_];
            unsigned long order = face_order_in_tet[0][bdy_];
            unsigned long face_tetra = nSide_*ik+order;

            // 临时变量存储排序后的 barx_ 和 phi_barx_ 的值
            double temp_barx[3][nBarx_];  
            double temp_phi_barx[nPhi_][nBarx_]; 

            // 读取原数组的值到临时变量，并将排序后的值写回原数组
            for (int j = 0; j < nBarx_; ++j) {  // 遍历每个点
                unsigned long original_index = indices[j];  // 排序前的位置索引
                for (int n = 0; n < 3; ++n) {  
                    temp_barx[n][j] = barx_[n][nBarx_ * bdy_ + original_index];
                }
                for (int n = 0; n < nPhi_; ++n) { 
                    temp_phi_barx[n][j] = phi_barx_[n][nBarx_ * face_tetra + original_index];
                }
            }
            // 交换ik面的积分点的顺序
            for (int j = 0; j < nBarx_; ++j) {  // 遍历每个点
                for (int n = 0; n < 3; ++n) { 
                    barx_[n][nBarx_ * bdy_ + j] = temp_barx[n][j];
                }
                for (int n = 0; n < nPhi_; ++n) {
                    phi_barx_[n][nBarx_ * face_tetra + j] = temp_phi_barx[n][j];
                }
            }
        }
    }

}


// 计算初值条件
void LinearDGSolver_3D_CycleBoundary::computeInitialCondition(unsigned long n, double *x, double *y, double *z, double **rtval) {
    double** u0 = new double*[nVars_];  
    for (int i = 0; i < nVars_; i++) {
        u0[i] = new double[dof_]{0.0};  // 初始化为 0
    }
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

    for (int i = 0; i < nVars_; i++) {
        delete[] u0[i];
    }
    delete[] u0;
}


// 计算单元边界数值流通量
void LinearDGSolver_3D_CycleBoundary::computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[3][5]) {   
    unsigned long** faceInfo = tetmesh_->face_info();
    unsigned long** tet_face_conn = tetmesh_->tet_face_connection();
    unsigned long** face_order_in_tet = tetmesh_->face_order_in_tet();

    // 边编号
    unsigned long idx_edge = tet_face_conn[ie][ik];

    // 外法向量
    double normal[3];
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
    double exactSolution[nVertex];
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
