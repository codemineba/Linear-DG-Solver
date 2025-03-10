#include "LinearDGSolver3D.h"


//计算基函数及其导数
void LinearDGSolver_3D::phi(unsigned long ik, double x, double y, double z, double *var, double *var_x, double *var_y, double *var_z){
    
    unsigned long **tet = tetmesh_->tetrahedron();
    double *x_=tetmesh_->x_coord();
    double *y_=tetmesh_->y_coord();
    double *z_=tetmesh_->z_coord();

    // 顶点坐标
    double x0 = x_[tet[0][ik]], x1 = x_[tet[1][ik]], x2 = x_[tet[2][ik]], x3 = x_[tet[3][ik]];
    double y0 = y_[tet[0][ik]], y1 = y_[tet[1][ik]], y2 = y_[tet[2][ik]], y3 = y_[tet[3][ik]];
    double z0 = z_[tet[0][ik]], z1 = z_[tet[1][ik]], z2 = z_[tet[2][ik]], z3 = z_[tet[3][ik]];
    double v0[dim_] = {x0, y0, z0};
    double v1[dim_] = {x1, y1, z1};
    double v2[dim_] = {x2, y2, z2};
    double v3[dim_] = {x3, y3, z3};
    double v[dim_]  = {x, y, z};

    // 体积坐标
    double volume=directed_3D_tetrahedron_volume(v0, v1, v2, v3);
    double volume0=directed_3D_tetrahedron_volume(v, v1, v2, v3);
    double volume1=directed_3D_tetrahedron_volume(v0, v, v2, v3);
    double volume2=directed_3D_tetrahedron_volume(v0, v1, v, v3);
    double eta0 = volume0/volume;
    double eta1 = volume1/volume;
    double eta2 = volume2/volume;

    var[0] = 1;
    var[1] = -1.0/4 + eta0;
    var[2] = -1.0/3 + eta0/3 + eta1;
    var[3] = -1.0/2 + eta0/2 + eta1/2 + eta2;

    if (var_x!=nullptr && var_y!=nullptr && var_z!=nullptr){
        double volume0_v0[3], volume1_v1[3], volume2_v2[3];
        directed_3D_tetrahedron_volume_gradient(v, v1, v2, v3, 0, volume0_v0);
        directed_3D_tetrahedron_volume_gradient(v0, v, v2, v3, 1, volume1_v1);
        directed_3D_tetrahedron_volume_gradient(v0, v1, v, v3, 2, volume2_v2);

        double tem_volumn_v[dim_][nPhi_];
        for(int i=0; i<dim_; i++){
            tem_volumn_v[i][0]=0.0;
            tem_volumn_v[i][1]=volume0_v0[i]/volume;
            tem_volumn_v[i][2]=((1.0/3.0)*volume0_v0[i]+volume1_v1[i])/volume;
            tem_volumn_v[i][3]=((1.0/2.0)*volume0_v0[i]+(1.0/2.0)*volume1_v1[i]+volume2_v2[i])/volume;
        }
        for(int i=0; i<nPhi_; i++){
            var_x[i]=tem_volumn_v[0][i];
            var_y[i]=tem_volumn_v[1][i];
            var_z[i]=tem_volumn_v[2][i];
        }
    }
}

// 计算通量, 返回3x5的矩阵f
void LinearDGSolver_3D::flux(double* u, double f[3][5]){

    double rho=u[0];   // 密度
    double v1=u[1]/u[0];  // x方向速度
    double v2=u[2]/u[0];  // y方向速度
    double v3=u[3]/u[0];  // z方向速度
    double v[dim_]={v1, v2, v3}; 
    double E=u[4];  // 能量
    double p =PRESSURE(E, rho, v); // 压强
    
    f[0][0] = rho*v1;        f[1][0] = rho*v2;        f[2][0] = rho*v3;
    f[0][1] = rho*v1*v1+p;   f[1][1] = rho*v1*v2;     f[2][1] = rho*v1*v3;
    f[0][2] = rho*v2*v1;     f[1][2] = rho*v2*v2+p;   f[2][2] = rho*v2*v3; 
    f[0][3] = rho*v3*v1;     f[1][3] = rho*v3*v2;     f[2][3] = rho*v3*v3+p;
    f[0][4] = (E+p)*v1;      f[1][4] = (E+p)*v2;      f[2][4] = (E+p)*v3;
    
}

void LinearDGSolver_3D::numerical_flux(double* u1, double *u2, double *normal, double* f){
    double fu1[dim_][5], fu2[dim_][5], fn[nVars_];
    flux(u1, fu1);
    flux(u2, fu2);
    for(int i=0; i<nVars_; i++){
        fn[i]=(fu1[0][i]+fu2[0][i])*normal[0]+(fu1[1][i]+fu2[1][i])*normal[1]+(fu1[2][i]+fu2[2][i])*normal[2];
    }

    // 计算谱半径
    double rho1=u1[0], E1=u1[4], v1[dim_]={u1[1]/u1[0], u1[2]/u1[0], u1[3]/u1[0]};
    double cs1 =SOUND_SPEED(PRESSURE(E1, rho1, v1), rho1);
    double lambda_u1 = abs(dot(v1, normal, dim_))+cs1;

    double rho2=u2[0], E2=u2[4], v2[dim_]={u2[1]/u2[0], u2[2]/u2[0], u2[3]/u2[0]};
    double cs2 =SOUND_SPEED(PRESSURE(E2, rho2, v2), rho2);    
    double lambda_u2 = abs(dot(v2, normal, dim_))+cs2;

    double alpha= lambda_u1 > lambda_u2 ? lambda_u1 : lambda_u2;


    for (int i=0; i<nVars_; i++){
        f[i]=(fn[i]-alpha*(u2[i]-u1[i]))/2.0;
    }
    
}


// 计算单元属性
void LinearDGSolver_3D::computeElementProperties(){
    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();
    unsigned long ** tet=tetmesh_->tetrahedron();
    unsigned long ** faceidx=tetmesh_->face_info();
    unsigned long ** tet_edge_conn=tetmesh_->tet_face_connection();

    // 计算体积
    for (unsigned long i=0; i<nElement; ++i){
        double v0[dim_]={x[tet[0][i]], y[tet[0][i]], z[tet[0][i]]};
        double v1[dim_]={x[tet[1][i]], y[tet[1][i]], z[tet[1][i]]};
        double v2[dim_]={x[tet[2][i]], y[tet[2][i]], z[tet[2][i]]};
        double v3[dim_]={x[tet[3][i]], y[tet[3][i]], z[tet[3][i]]};

        volume_[i]=abs(directed_3D_tetrahedron_volume(v0, v1, v2, v3));
    }

    // 计算单元内求积节点
    double alpha=0.58541020, beta=0.13819660;
    double hateta[nHatx_][4]={  // 标准单元内的4个求积节点
        {alpha, beta, beta, 1-alpha-beta-beta},
        {beta, alpha, beta, 1-beta-alpha-beta},
        {beta, beta, alpha, 1-beta-beta-alpha},
        {beta, beta, beta,  1-beta-beta-beta}
    };
    for(unsigned long i=0; i<nElement; ++i){ //按单元循环
        for(int j=0; j<nHatx_; ++j){ // 按四个求积点循环
        hatx_[0][nHatx_*i+j]=hatx_[1][nHatx_*i+j]=hatx_[2][nHatx_*i+j]=0.0;
            for(int k=0; k<4; ++k){  // 按四个顶点循环
                hatx_[0][nHatx_*i+j]+=x[tet[k][i]]*hateta[j][k];
                hatx_[1][nHatx_*i+j]+=y[tet[k][i]]*hateta[j][k];
                hatx_[2][nHatx_*i+j]+=z[tet[k][i]]*hateta[j][k];
            }
        }
    }

    double bareta[nBarx_][3]={
        {2.0/3.0, 1.0/6.0, 1.0/6.0},
        {1.0/6.0, 2.0/3.0, 1.0/6.0},
        {1.0/6.0, 1.0/6.0, 2.0/3.0}
    };
    for (unsigned long i=0; i<nFace; ++i){  // 按面循环
        for(int j=0; j<nBarx_; j++){  // 按三个求积点循环
        barx_[0][nBarx_*i+j]=barx_[1][nBarx_*i+j]=barx_[2][nBarx_*i+j]=0.0;
            for(int k=0; k<3; k++){  // 按三角形顶点循环
                barx_[0][nBarx_*i+j]+=x[faceidx[k][i]]*bareta[j][k];
                barx_[1][nBarx_*i+j]+=y[faceidx[k][i]]*bareta[j][k];
                barx_[2][nBarx_*i+j]+=z[faceidx[k][i]]*bareta[j][k];
            }
        }
        double v0[3]={x[faceidx[0][i]], y[faceidx[0][i]], z[faceidx[0][i]]}, v1[3]={x[faceidx[1][i]], y[faceidx[1][i]], z[faceidx[1][i]]}, v2[3]={x[faceidx[2][i]], y[faceidx[2][i]], z[faceidx[2][i]]};
        face_area_[i]=abs(directed_3D_triangle_area(v0, v1, v2));
    }

    // 计算表面积
    for (unsigned long i=0; i<nElement; ++i){
        sufacearea_[i] = 0.0;
        for(int j=0; j<nSide_; j++){  // 4个面
        sufacearea_[i] += face_area_[tet_edge_conn[j][i]];
        }
    }

}


// 计算外法向量
void LinearDGSolver_3D::computeOuterNormal(){
    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();
    unsigned long ** tet=tetmesh_->tetrahedron();
    unsigned long** tet_edge_conn = tetmesh_->tet_face_connection();
    unsigned long** faceInfo = tetmesh_->face_info();

    for(unsigned long ik=0; ik<nElement; ik++){
        for(int ie=0; ie<nSide_; ie++){  // 每个面 
            unsigned long p0, p1, p2, p3;

            p1 = faceInfo[0][tet_edge_conn[ie][ik]];
            p2 = faceInfo[1][tet_edge_conn[ie][ik]];
            p3 = faceInfo[2][tet_edge_conn[ie][ik]]; 
            // p0 是 tet[0][ik], tet[1][ik], tet[2][ik], tet[3][ik] 中 p1, p2, p3 除外的
            p0 = (tet[0][ik] != p1 && tet[0][ik] != p2 && tet[0][ik] != p3) ? tet[0][ik] :
                   (tet[1][ik] != p1 && tet[1][ik] != p2 && tet[1][ik] != p3) ? tet[1][ik] :
                   (tet[2][ik] != p1 && tet[2][ik] != p2 && tet[2][ik] != p3) ? tet[2][ik] :
                   tet[3][ik];

            // 计算法向量
            double a1a2[dim_]={x[p2]-x[p1],y[p2]-y[p1],z[p2]-z[p1]};
            double a1a3[dim_]={x[p3]-x[p1],y[p3]-y[p1],z[p3]-z[p1]};
            double a1a0[dim_]={x[p0]-x[p1],y[p0]-y[p1],z[p0]-z[p1]};
            double a1a2_cross_a1a3[dim_];
            cross(a1a2, a1a3, a1a2_cross_a1a3);
            double magnitude=l2_norm(a1a2_cross_a1a3, dim_);
            if(dot(a1a2_cross_a1a3, a1a0, dim_)<0){
                outerNormal_[0][ik*nSide_+ie]=a1a2_cross_a1a3[0]/magnitude;
                outerNormal_[1][ik*nSide_+ie]=a1a2_cross_a1a3[1]/magnitude;
                outerNormal_[2][ik*nSide_+ie]=a1a2_cross_a1a3[2]/magnitude;
            }else if(dot(a1a2_cross_a1a3, a1a0, dim_)>0){
                outerNormal_[0][ik*nSide_+ie]=-a1a2_cross_a1a3[0]/magnitude;
                outerNormal_[1][ik*nSide_+ie]=-a1a2_cross_a1a3[1]/magnitude;
                outerNormal_[2][ik*nSide_+ie]=-a1a2_cross_a1a3[2]/magnitude;
            }
        }
    }
}


// 计算基函数在相关节点上的值
void LinearDGSolver_3D::computeBasisOnNodes(){
    unsigned long** tet_edge_conn = tetmesh_->tet_face_connection();

    // 计算基函数相关值 
    // 在不同单元顶点基函数值是相同的
    phi_vec_[0][0] = 1.0;     phi_vec_[0][1] = 1.0;     phi_vec_[0][2] = 1.0;     phi_vec_[0][3] = 1.0; 
    phi_vec_[1][0] = 3.0/4.0; phi_vec_[1][1] =-1.0/4.0; phi_vec_[1][2] =-1.0/4.0; phi_vec_[1][3] =-1.0/4.0;
    phi_vec_[2][0] = 0.0;     phi_vec_[2][1] = 2.0/3.0; phi_vec_[2][2] =-1.0/3.0; phi_vec_[2][3] =-1.0/3.0;
    phi_vec_[3][0] = 0.0;     phi_vec_[3][1] = 0.0;     phi_vec_[3][2] = 1.0/2.0; phi_vec_[3][3] =-1.0/2.0;
    for(unsigned long ik=0; ik<nElement; ik++){ //按单元循环
        for(int ie=0; ie<nSide_; ie++){ //按面循环
            double phival[nPhi_], phival_x[nPhi_], phival_y[nPhi_], phival_z[nPhi_];
            unsigned long faceIndex = tet_edge_conn[ie][ik];
            phi(ik, hatx_[0][ik*nHatx_+ie], hatx_[1][ik*nHatx_+ie], hatx_[2][ik*nHatx_+ie], phival, phival_x, phival_y, phival_z);
            for(int j=0; j<nPhi_; j++){
                phi_hatx_[j][ik*nHatx_+ie] = phival[j];    
                phi_hatx_x_[j][ik*nHatx_+ie] = phival_x[j]; 
                phi_hatx_y_[j][ik*nHatx_+ie] = phival_y[j];
                phi_hatx_z_[j][ik*nHatx_+ie] = phival_z[j];   // nHatx_ = nSide_ 
            }
            for(int m=0; m<nBarx_; m++){
                phi(ik, barx_[0][nBarx_*faceIndex+m], barx_[1][nBarx_*faceIndex+m], barx_[2][nBarx_*faceIndex+m], phival);
                for(int j=0; j<nPhi_; j++){
                    phi_barx_[j][nBarx_*(ik*nSide_+ie)+m] = phival[j]; 
                }
            }
        }
    }
}


// 计算初始状态
void LinearDGSolver_3D::computeInitialData() {
    for (int i=0; i<nVars_; i++){
        u0_[i] = new double[dof_];
    }

    // 初值条件
    computeInitialCondition(dof_, hatx_[0], hatx_[1], hatx_[2], u0_);

    // 赋值给u_
    for (unsigned long i = 0; i<nElement ; i++){
        for (int j = 0; j<nPhi_; j++) {
            for (int k=0; k<nVars_; k++){
                u_[k][i*nPhi_+j]=u0_[k][i*nPhi_+j];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
            }
        }
    }
}


void LinearDGSolver_3D::computeFluxOnElement(unsigned long ik, double** u, double f[3][20]){
    double fu[dim_][5]; // 用于储存该边的f(u)通量
    double ui[nVars_];  // 储存该边的物理量
    for(int i=0; i<nHatx_; i++){  // 按四个求积点循环
        for(int n=0; n<nVars_; n++){  // 5个物理量
            ui[n]=0.0;
            for(int j=0; j<nPhi_; j++){  // 四个基函数
                ui[n] += u[n][ik*nPhi_+j]*phi_hatx_[j][ik*nHatx_+i];  // 计算数值解在求积点上的值
            }
        }
        flux(ui, fu);  // 计算该边通
        // 将该边通量储存到f中 f即储存了一个单元四个通量
        for(int k=0; k<dim_; k++){
            for(int n=0; n<nVars_; n++){
                f[k][i*nVars_+n] = fu[k][n];
            }
        }
    }
}


void LinearDGSolver_3D::computeNumericalSolutionOnGaussPoints(unsigned long ik, unsigned long ie, double** u, double uh[3][5]){
    for (int i=0; i<nVars_; i++){      // 遍历5个物理量
        for (int m=0; m<nBarx_; m++){
            uh[m][i]=0;
            for (int j=0; j<nPhi_; j++){  // 遍历基函数
                uh[m][i] += u[i][ik*nPhi_+j] * phi_barx_[j][nBarx_*(ik*nSide_+ie)+m];
            }
        }
    }
}



// 计算单元内部流通量, 返回2x4矩阵
void LinearDGSolver_3D::computeNumericalFluxOnInteriorEdge(unsigned long ik, unsigned long ie, double** u, double f[3][5]){
    
    unsigned long** faceInfo = tetmesh_->face_info();
    unsigned long** tet_face_conn = tetmesh_->tet_face_connection();
    unsigned long** face_order_in_tet = tetmesh_->face_order_in_tet();

    // 边编号
    unsigned long idx_face = tet_face_conn[ie][ik];    

    // 外法向量
    double normal[dim_];
    normal[0] = outerNormal_[0][ik*nSide_+ie];
    normal[1] = outerNormal_[1][ik*nSide_+ie];
    normal[2] = outerNormal_[2][ik*nSide_+ie];

    // 该面相邻的两个四面体
    unsigned long k1 = faceInfo[3][idx_face];
    unsigned long k2 = faceInfo[4][idx_face];

    // idx_face对应的另一个三角形ik_ 以及 idx_face在ik_中的order ie_
    unsigned long ik_, ie_; 
    (k1 == ik) ? (ik_ = k2, ie_ = face_order_in_tet[1][idx_face]) : 
    (k2 == ik) ? (ik_ = k1, ie_ = face_order_in_tet[0][idx_face]) :
    (ik_ = nElement+1, ie_ = 4);
    if (ik_ == nElement+1) {  // 找不到 (边界)
        std::cerr << "Error: can't find another tetrahedron" << std::endl;
        throw -1;
    }

    // 计算u(x, y) 在边界上积分点的值
    double u1[nBarx_][5], u2[nBarx_][5];
    computeNumericalSolutionOnGaussPoints(ik, ie, u, u1);
    computeNumericalSolutionOnGaussPoints(ik_, ie_, u, u2);

    // 计算流通量
    for (int m=0; m<nBarx_; m++){  // 对于四个积分点
        double ff[nVars_];
        numerical_flux(u1[m], u2[m], normal, ff);
        for(int n=0; n<nVars_; n++){ // 对于五个物理量
            f[m][n] = ff[n]; 
        }
    }
}
    



// 计算整个单元边界上的数值流通量, 返回6x4矩阵
void LinearDGSolver_3D::computeNumericalFluxOnElement(unsigned long ik, double** u, double f[12][5]){

    unsigned long** faceInfo = tetmesh_->face_info();
    unsigned long **tet_face_conn_ = tetmesh_->tet_face_connection();
    
    double flux[dim_][5];
    for(int i=0; i<nSide_; i++){  // 四个面
        unsigned long ie = tet_face_conn_[i][ik];  // 该面索引
        // 判断该边是否是边界
        if (faceInfo[4][ie] == nElement+1){  // 是边界
            computeNumericalFluxOnBoundary(ik, i, u, flux);
        } 
        else{
            computeNumericalFluxOnInteriorEdge(ik, i, u, flux);
        }
        
        // 将该边通量储存到f中 f即储存了一个单元四个面(12个点)的数值通量
        for(int k=0; k<nBarx_; k++){
            for(int j=0; j<nVars_; j++){
                f[nBarx_*i+k][j] = flux[k][j];
            }
        }
    }
}


// 计算空间离散, 返回 4x dof_ 的矩阵f
void LinearDGSolver_3D::computeSpaceDiscretization(double** u, double**f) {

    unsigned long** conn=tetmesh_->tet_face_connection();
    for(int n=0; n<nVars_; n++){
        for (unsigned long j=0; j<dof_; j++){
            f[n][j]=0.0;
        }
    }
    double fu[dim_][20], flux[nBarx_*nSide_][5];  // 20=nVars_*nHatx_  5=nVars_
    for (unsigned long ik=0; ik<nElement; ik++){  // 单元数
        computeFluxOnElement(ik, u, fu);  // 计算k单元的f(u)通量
        computeNumericalFluxOnElement(ik, u, flux);  // 计算k单元的数值通量
        for (int j=0; j<nPhi_; j++){  // 基函数循环
            double phiProduct=0.0;  // 基函数内积
            if(j==0){phiProduct=1.0/6.0;}
            else if(j==1){phiProduct=1.0/160.0;}
            else if(j==2){phiProduct=1.0/180.0;}
            else if(j==3){phiProduct=1.0/240.0;}
            double aj = 6.0*volume_[ik]*phiProduct;

            // 常微分方程的第一项
            double scale1 = volume_[ik]/(nHatx_*aj); 
            for (int i=0; i<nHatx_; i++){ // 积分点循环
                for(int n=0; n<nVars_; n++){   // 对于5个物理量
                    f[n][ik*nPhi_+j] += scale1 * (fu[0][i*nVars_+n] * phi_hatx_x_[j][ik*nHatx_+i] + fu[1][i*nVars_+n] * phi_hatx_y_[j][ik*nHatx_+i]
                    + fu[2][i*nVars_+n] * phi_hatx_z_[j][ik*nHatx_+i]); 
                }   
            }
            // 常微分方程的第二项
            for (int ie=0; ie<nSide_; ie++){   // 边循环
                double scale2 = face_area_[conn[ie][ik]]/(nBarx_*aj);
                for (int m=0; m<nBarx_; m++){  // 面上求积点循环
                    for(int n=0; n<nVars_; n++){  // 对于5个物理量
                        f[n][ik*nPhi_+j] -= scale2 * flux[nBarx_*ie+m][n] * phi_barx_[j][nBarx_*(ik*nSide_+ie)+m];
                    }
                }
            }
        }
    }
    
}

double LinearDGSolver_3D::computeTimeStep(double** u){
    
    double dm[nElement]={0}; 

    for(unsigned long i=0; i<nElement; i++){
        double u_average[nVars_];
        for(int j=0; j<nVars_; j++){
            u_average[j] = (u[j][i*nPhi_] + u[j][i*nPhi_+1] + u[j][i*nPhi_+2] + u[j][i*nPhi_+3]) / nPhi_;
        }
        double rho=u_average[0], E=u_average[4], v[dim_]={u_average[1]/u_average[0], u_average[2]/u_average[0], u_average[3]/u_average[0]};
        double cs =SOUND_SPEED(PRESSURE(E, rho, v), rho);
        dm[i] = (l2_norm(v, dim_) + cs) * (sufacearea_[i] / volume_[i]);
    }

    int n = sizeof(dm) / sizeof(dm[0]); // 数组大小 (nElement)
    double max_dm = *std::max_element(dm, dm + n);  // 寻找最大值 

    return CFL/max_dm;
    
}


void LinearDGSolver_3D::computeTimeDiscretization(double total_time, double* save_time_points, int save_time_num, double dt){
    
    // dt默认为0 不提供参数 即采取变时间步长
    if (dt)
        timeInterval_ = dt;

    // save_time_num 默认唯一 即默认只保存第一步的值和最后一步的值
    for (int i=0; i<nVars_; i++){
        u_[i]=new double[dof_*save_time_num];
    }

    // 计算初始状态
    computeInitialData();

    double *u_now[nVars_]; // 储存当前的值
    double *u_middle[nVars_]; // 储存中间变量
    double *L_jk[nVars_]; // 储存用于空间离散的值
    for (int i=0; i<nVars_; i++){  // 初始化数组
        u_now[i]=new double[dof_];
        u_middle[i]=new double[dof_];
        L_jk[i]=new double[dof_];
    }
    for (int m=0; m<nVars_; m++){   // 每个物理量
        for (unsigned long i=0; i<dof_; i++) {  // 长度为dof_
            u_now[m][i] = u_[m][i];
        }
    }

    SynchronizationUpdate(u_now);  // 同步更新    

    // 输出日志信息
    std::cout << "time_step " << 0 << ", " << LogMessage() << std::endl;

    int time_step_num=1, time_step_index=1;   
    num_times_point_save_ = save_time_num;
    times_point_save_ = new double[num_times_point_save_];
    if(save_time_points){
        for(int i = 0; i<num_times_point_save_; i++){
            times_point_save_[i] = save_time_points[i];
        }
        // 对要保存的时间点进行排序,方便后面记录
        std::sort(save_time_points, save_time_points + num_times_point_save_); 
        time_step_index = (int(save_time_points[0]) == 0)? 1 : 0;  // 如果选择保存初始值, 就从索引1开始保存 (索引0留给初始值)
    }else{ 
        times_point_save_[0] = 0;   // 默认只保留开始值和终止值
        times_point_save_[1] = total_time;
    }

    while(1){  // 但计算实时时间达到(接近)total_time即跳出循环

        // 更新下一阶段的时间步长
        if (!dt)
        timeInterval_ = computeTimeStep(u_now);

        // 二阶 TVD Runng-kutta 计算格式
        computeSpaceDiscretization(u_now, L_jk);
        // 计算中间值 第（n + 1/2）步
        for (int m=0; m<nVars_; m++){   // 每个物理量
            for (unsigned long i=0; i<dof_; i++){
                u_middle[m][i] = u_now[m][i] + timeInterval_ * L_jk[m][i]; 
            }
        }
        
        computeSpaceDiscretization(u_middle, L_jk);
        // 计算 第n+1步
        for (int m=0; m<nVars_; m++){   // 每个物理量
            for (unsigned long i=0; i<dof_; i++){
                u_now[m][i] = 0.5 * (u_now[m][i] + u_middle[m][i]) + 0.5 * timeInterval_ * L_jk[m][i]; 
            }
        }

        timePoint_ += timeInterval_;  //更新实时时刻
        SynchronizationUpdate(u_now); // 同步更新

        // 输出日志信息
        std::cout << "time_step " << time_step_num << ", " << LogMessage() << std::endl;
        
        // 检查是否需要保存该时间点的结果
        if(save_time_points && time_step_index<num_times_point_save_){
            double save_time_point = save_time_points[time_step_index];
            double interval_save = abs(save_time_point - timePoint_);  // 当前与目标保存时长的时间距离
            double intervalNext_save = abs(save_time_point - (timePoint_ + timeInterval_)); // 下一阶段与目标保存时长的时间距离
            int shift = dof_ * time_step_index;
            if (interval_save < intervalNext_save) {  // 找出与目标保存时长距离最近的时间点
                std::cout<<"Saved, time is " << timePoint_ << std::endl;
                for (int m=0; m<nVars_; m++){   // 每个物理量
                    for (unsigned long i=0; i<dof_; i++) {  // 长度为dof_
                        u_[m][i+shift] = u_now[m][i];
                    }
                }
                time_step_index ++;
            }  
        }

        // 检查是否已到达目标计算时长
        double interval_total = abs(total_time - timePoint_);  // 当前与目标计算时长的时间距离
        double intervalNext_total = abs(total_time - (timePoint_ + timeInterval_)); // 下一阶段与目标计算时长的时间距离
        if (interval_total < intervalNext_total) {  // 找出与目标计算时长距离最近的时间点
            std::cout<<"Iterated over! Total time is " << timePoint_ << std::endl;
            if(!save_time_points){  // 默认保存最后一步的值
                for (int m=0; m<nVars_; m++){   // 每个物理量
                    for (unsigned long i=0; i<dof_; i++) {  // 长度为dof_
                        u_[m][i+dof_] = u_now[m][i];
                    }
                }
            }
            break;
        }   
        time_step_num ++;  // 否则, 继续时间离散
    }
     // 释放空间
    for (int i = 0; i < nVars_; i++) {
        delete[] u_now[i];
        delete[] u_middle[i];
        delete[] L_jk[i];
    }
}

// 计算某个时刻的单元顶点的一个数值解 
void LinearDGSolver_3D::NumericalSolutionOnVertex(double** u, int idx, double* var){
    unsigned long **conn=tetmesh_->tetrahedron();
    
    unsigned int count[nVertex];
    for (unsigned long j=0; j<nVertex; j++){  // 更新
        count[j]=0;
        var[j]=0.0;
    }
    for (unsigned long ik=0; ik<nElement; ik++){
        for (int i=0; i<nSide_; i++){  // 按顶点循环
            for (int j=0; j<nPhi_; j++){   // 基函数循环
                var[conn[i][ik]]+= u[idx][ik*nPhi_+j]*phi_vec_[j][i];   
            }
            count[conn[i][ik]]++;
        }
    }
    for (unsigned long j=0; j<nVertex; j++){
        var[j]=var[j]/count[j];
    }
}

// 计算某个时刻的单元顶点的全部数值解
void LinearDGSolver_3D::NumericalSolutionOnVertex(double** u, double** var){

    unsigned long **conn=tetmesh_->tetrahedron();

    unsigned int count[nVertex];
    for (unsigned long j=0; j<nVertex; j++){
        count[j]=0;
        for(int k=0; k<nVars_; k++){
            var[k][j] =0.0;
        }
    }
    for (unsigned long ik=0; ik<nElement; ik++){  // 按单元循环
        for (int i=0; i<nSide_; i++){  // 按顶点循环
            for (int j=0; j<nPhi_; j++){   // 基函数循环
                for (int n=0; n<nVars_; n++){
                    var[n][conn[i][ik]]+= u[n][ik*nPhi_+j]*phi_vec_[j][i];   
                }
            }
            count[conn[i][ik]]++;
        }
    }
    for (int i=0; i<nVars_; i++){
        for (unsigned long j=0; j<nVertex; j++){
            var[i][j]=var[i][j]/count[j];
        }
    }
}



// 用于画图
void LinearDGSolver_3D::SynchronizationUpdate(double** u){
    NumericalSolutionOnVertex(u, 0, rho);
}

    

void LinearDGSolver_3D::outputTecPlotDataFile(const std::string &fname, double time_point){

    double *x=tetmesh_->x_coord();
    double *y=tetmesh_->y_coord();
    double *z=tetmesh_->z_coord();
    unsigned long **tet=tetmesh_->tetrahedron();

    //寻找记录的保存时间点
    int index = -1;
    for(int i =0; i<num_times_point_save_; i++){
        if(times_point_save_[i] == time_point){
            index = i;
            break;
        }
    }
    // 如果未找到匹配的时间点，输出错误提示并终止函数
    if (index == -1) {
        std::cerr << "Error: Time_point " << time_point 
             << " not found in saved_time_points_" << std::endl;
        throw -1;
    }

    int shift = index*dof_;
    double* uj[nVars_];
    for (int i=0; i<nVars_; i++){
        uj[i]=new double[dof_];
        for (unsigned long j=0; j<dof_; j++){
            uj[i][j]=u_[i][shift + j];
        }
    }

    // 用于保留单元顶点上的数值解
    double* uh[nVars_];
    for (int i=0; i<nVars_; i++){
        uh[i]=new double[nVertex];
        for (unsigned long j=0; j<nVertex; j++){
            uh[i][j]=0.0;
        }
    }
    
    // 计算time_step时的单元顶点上的数值解
    NumericalSolutionOnVertex(uj, uh);

    std::ofstream output(fname);
    output << "TITLE=\"Tetrahedral mesh\"" << std::endl;
    output << "VARIABLES= \"x\" \"y\" \"z\" \"rho\" \"vx\" \"vy\" \"vz\" \"E\" \"p\"" << std::endl;
    output << "ZONE T=\"none\","
           << "N=" << nVertex << ","
           << "E=" << nElement << ","
           << "ET=TETRAHEDRON, F=FEPOINT" << std::endl;
    
    for (unsigned long i = 0; i < nVertex; i++) {
        double rho = uh[0][i];
        double vx = uh[1][i] / rho;
        double vy = uh[2][i] / rho;
        double vz = uh[3][i] / rho;
        double v[dim_] = {vx, vy, vz};
        double E = uh[4][i];
        double p = PRESSURE(E, rho, v);
    
        output << std::setprecision(15) << std::setw(20) << x[i] << "\t";
        output << std::setprecision(15) << std::setw(20) << y[i] << "\t";
        output << std::setprecision(15) << std::setw(20) << z[i] << "\t";
        output << std::setprecision(15) << std::setw(20) << rho << "\t";
        output << std::setprecision(15) << std::setw(20) << vx << "\t";
        output << std::setprecision(15) << std::setw(20) << vy << "\t";
        output << std::setprecision(15) << std::setw(20) << vz << "\t";
        output << std::setprecision(15) << std::setw(20) << E << "\t";
        output << std::setprecision(15) << std::setw(20) << p << std::endl;
    }
    
    for (unsigned long i = 0; i < nElement; i++) {
        output << std::setw(10) << tet[0][i] + 1 << "\t";
        output << std::setw(10) << tet[1][i] + 1 << "\t";
        output << std::setw(10) << tet[2][i] + 1 << "\t";
        output << std::setw(10) << tet[3][i] + 1 << std::endl;
    }
    
    output.close();

    std::cout << "successfully written to " << fname << std::endl;

    for (int i = 0; i < nVars_; i++){
        delete[] uh[i];
        delete[] uj[i];
    }
}

void LinearDGSolver_3D::outputVTKDataFile(const std::string &fname, double time_point) {
    double *x = tetmesh_->x_coord();
    double *y = tetmesh_->y_coord();
    unsigned long **conn = tetmesh_->tetrahedron();

    //寻找记录的保存时间点
    int index = -1;
    for(int i =0; i<num_times_point_save_; i++){
        if(times_point_save_[i] == time_point){
            index = i;
            break;
        }
    }
    // 如果未找到匹配的时间点，输出错误提示并终止函数
    if (index == -1) {
        std::cerr << "Error: Time_point " << time_point 
             << " not found in saved_time_points_" << std::endl;
        throw -1;
    }

    int shift = index*dof_;
    double* uj[nVars_];
    for (int i=0; i<nVars_; i++){
        uj[i]=new double[dof_];
        for (unsigned long j=0; j<dof_; j++){
            uj[i][j]=u_[i][shift + j];
        }
    }

    // 用于保留单元顶点上的数值解
    double *uh[nVars_];
    for (int i = 0; i < nVars_; i++) {
        uh[i] = new double[nVertex];
        for (unsigned long j = 0; j < nVertex; j++) {
            uh[i][j] = 0.0;
        }
    }

    // 计算 time_step 时的单元顶点上的数值解
    NumericalSolutionOnVertex(uj, uh);

    std::ofstream output(fname);
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "Triangular mesh" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // 写入点坐标
    output << "POINTS " << nVertex << " double" << std::endl;
    for (unsigned long i = 0; i < nVertex; i++) {
        output << std::setprecision(15) << std::setw(20) << x[i] << " ";
        output << std::setprecision(15) << std::setw(20) << y[i] << " 0.0" << std::endl; // z 方向设为 0.0
    }

    // 写入单元连接信息（三角形）
    output << "CELLS " << nElement << " " << nElement * 4 << std::endl; // 每个三角形有 3 个顶点，行开头有个计数 3
    for (unsigned long i = 0; i < nElement; i++) {
        output << "3 " << conn[0][i] << " " << conn[1][i] << " " << conn[2][i] << std::endl;
    }

    // 写入单元类型（三角形类型在 VTK 中为 5）
    output << "CELL_TYPES " << nElement << std::endl;
    for (unsigned long i = 0; i < nElement; i++) {
        output << "5" << std::endl; // VTK_TRIANGLE 类型对应编号为 5
    }

    // 写入顶点上的标量数据 rho, vx, vy, E 和 p
    output << "POINT_DATA " << nVertex << std::endl;

    // 写入 rho 数据
    output << "SCALARS rho double 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (unsigned long i = 0; i < nVertex; i++) {
        double rho = uh[0][i];
        output << std::setprecision(15) << std::setw(20) << rho << std::endl; // rho
    }

    // 写入速度向量数据 (vx, vy)
    output << "VECTORS velocity double" << std::endl;
    for (unsigned long i = 0; i < nVertex; i++) {
        double rho = uh[0][i];
        double vx = uh[1][i] / rho; // vx
        double vy = uh[2][i] / rho; // vy
        output << std::setprecision(15) << std::setw(20) << vx << " ";
        output << std::setprecision(15) << std::setw(20) << vy << " ";
        output << "0.0" << std::endl; // z 方向速度为 0.0
    }

    // 写入能量数据 E
    output << "SCALARS E double 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (unsigned long i = 0; i < nVertex; i++) {
        double E = uh[3][i];
        output << std::setprecision(15) << std::setw(20) << E << std::endl; // E
    }

    // 写入压力数据 p
    output << "SCALARS p double 1" << std::endl;
    output << "LOOKUP_TABLE default" << std::endl;
    for (unsigned long i = 0; i < nVertex; i++) {
        double rho = uh[0][i];
        double vx = uh[1][i] / rho;
        double vy = uh[2][i] / rho;
        double vz = uh[3][i] / rho;
        double v[dim_] = {vx, vy, vz};
        double E = uh[4][i];
        double p = PRESSURE(E, rho, v); // p 的计算
        output << std::setprecision(15) << std::setw(20) << p << std::endl;
    }

    output.close();

    // 释放内存
    for (int i = 0; i < nVars_; i++){
        delete[] uh[i];
        delete[] uj[i];
    }
}
