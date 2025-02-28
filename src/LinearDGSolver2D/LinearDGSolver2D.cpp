#include "LinearDGSolver2D.h"


//计算基函数及其导数
void LinearDGSolver_2D::phi(unsigned long ik, double x, double y, double *var, double *var_x, double *var_y){

    unsigned long **tri = trimesh_->triangle();
    double *x_=trimesh_->x_coord();
    double *y_=trimesh_->y_coord();
        // 顶点坐标
    double x0 = x_[tri[0][ik]];
    double x1 = x_[tri[1][ik]];
    double x2 = x_[tri[2][ik]];
    double y0 = y_[tri[0][ik]];
    double y1 = y_[tri[1][ik]];  
    double y2 = y_[tri[2][ik]];
    double v0[dim_] = {x0, y0};
    double v1[dim_] = {x1, y1};
    double v2[dim_] = {x2, y2};
    double v[dim_]  = {x, y};

    // 面积坐标
    double s=directed_2D_triangle_area(v0, v1, v2);
    double s0=directed_2D_triangle_area(v, v1, v2);
    double s1=directed_2D_triangle_area(v0, v, v2);
    double xi = s0/s;
    double eta = s1/s;

    var[0] = 1-2*xi;
    var[1] = 1-2*eta;
    var[2] = 1 - 2 * (1 - xi - eta);

    if (var_x!=nullptr && var_y!=nullptr){
        double s0_v0[dim_], s1_v1[dim_];   //关于v0的s0梯度 和 关于v1的s1梯度
        directed_2D_triangle_area_gradient(v, v1, v2, 0, s0_v0);
        directed_2D_triangle_area_gradient(v0, v, v2, 1, s1_v1);

        double tem_s_v[dim_][nPhi_];
        for(int i=0; i<dim_; i++){
            tem_s_v[i][0]=-2*s0_v0[i]/s;
            tem_s_v[i][1]=-2*s1_v1[i]/s;
            tem_s_v[i][2]=(2*s0_v0[i]+2*s1_v1[i])/s;
        }
        for(int i=0; i<nPhi_; i++){
            var_x[i]=tem_s_v[0][i];
            var_y[i]=tem_s_v[1][i];
        }
    }

    // unsigned long **tri_edge_connection = trimesh_->tri_edge_connection();

    // // 中点坐标
    // double hatx0 = hatx_[0][tri_edge_connection[1][ik]];
    // double hatx1 = hatx_[0][tri_edge_connection[2][ik]];
    // double hatx2 = hatx_[0][tri_edge_connection[0][ik]];
    // double haty0 = hatx_[1][tri_edge_connection[1][ik]];
    // double haty1 = hatx_[1][tri_edge_connection[2][ik]];  
    // double haty2 = hatx_[1][tri_edge_connection[0][ik]];

    // double hatv0[2] = {hatx0, haty0};
    // double hatv1[2] = {hatx1, haty1};
    // double hatv2[2] = {hatx2, haty2};
    // double v[2] = {x, y};

    // double hatS  = directed_2D_triangle_area(hatv0, hatv1, hatv2);
    // double hatS0 = directed_2D_triangle_area(v, hatv1, hatv2);
    // double hatS1 = directed_2D_triangle_area(hatv0, v, hatv2);
    // double hatS2 = directed_2D_triangle_area(hatv0, hatv1, v);

    // var[0] = hatS0/ hatS;
    // var[1] = hatS1/ hatS;
    // var[2] = hatS2/ hatS;

    // if (var_x!=nullptr){
    //     var_x[0] = (haty1 - haty2) / (hatS*2);
    //     var_x[1] = (haty2 - haty0) / (hatS*2);
    //     var_x[2] = (haty0 - haty1) / (hatS*2);

    //     var_y[0] = (hatx2 - hatx1) / (hatS*2);
    //     var_y[1] = (hatx0 - hatx2) / (hatS*2);
    //     var_y[2] = (hatx1 - hatx0) / (hatS*2);
    // }
}

// 计算通量, 返回2x4的矩阵f
void LinearDGSolver_2D::flux(double* u, double f[2][4]){

    double rho=u[0];   // 密度
    double v1=u[1]/u[0];  // x方向速度
    double v2=u[2]/u[0];  // y方向速度
    double v[dim_]={v1, v2}; 
    double E=u[3];  // 能量
    double p =PRESSURE(E, rho, v); // 压强
    
    f[0][0] = rho*v1,        f[1][0] = rho*v2,
    f[0][1] = rho*v1*v1+p,   f[1][1] = rho*v1*v2,
    f[0][2] = rho*v1*v2,     f[1][2] = rho*v2*v2+p,
    f[0][3] = (E+p)*v1,      f[1][3] = (E+p)*v2;
    
}

void LinearDGSolver_2D::numerical_flux(double* u1, double *u2, double *normal, double* f){
    double fu1[dim_][4], fu2[dim_][4], fn[nVars_];  // 第二维必须指定数字才能传入函数  nVars_=4 
    flux(u1, fu1);
    flux(u2, fu2);
    for(int i=0; i<nVars_; i++){
        fn[i]=(fu1[0][i]+fu2[0][i])*normal[0]+(fu1[1][i]+fu2[1][i])*normal[1];
    }

    // 计算谱半径
    double rho1=u1[0], E1=u1[3], v1[dim_]={u1[1]/u1[0], u1[2]/u1[0]};
    double cs1 =SOUND_SPEED(PRESSURE(E1, rho1, v1), rho1);
    double lambda_u1 = abs(dot(v1, normal, dim_))+cs1;

    double rho2=u2[0], E2=u2[3], v2[dim_]={u2[1]/u2[0], u2[2]/u2[0]};
    double cs2 =SOUND_SPEED(PRESSURE(E2, rho2, v2), rho2);
    double lambda_u2 = abs(dot(v2, normal, dim_))+cs2;

    double alpha= lambda_u1 > lambda_u2 ? lambda_u1 : lambda_u2;
    
    for (int i=0; i<nVars_; i++){
        f[i]=(fn[i]-alpha*(u2[i]-u1[i]))/2.0;
    }
}


// 计算单元属性
void LinearDGSolver_2D::computeElementProperties(){
    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    unsigned long ** tri=trimesh_->triangle();
    unsigned long ** edgeidx=trimesh_->edge_info();
    unsigned long ** tri_edge_conn=trimesh_->tri_edge_connection();

    // 计算面积
    for (unsigned long i=0; i<nElement; ++i){

        double v0[dim_]={x[tri[0][i]], y[tri[0][i]]};
        double v1[dim_]={x[tri[1][i]], y[tri[1][i]]};
        double v2[dim_]={x[tri[2][i]], y[tri[2][i]]};

        area_[i]=abs(directed_2D_triangle_area(v0, v1, v2));
    }

    // 计算边界边长和中点高斯点坐标
    double c=2.0*sqrt(3.0);
    for (unsigned long i=0; i<nEdge; ++i){
        hatx_[0][i]=(x[edgeidx[0][i]]+x[edgeidx[1][i]])/2;
        hatx_[1][i]=(y[edgeidx[0][i]]+y[edgeidx[1][i]])/2;
        double dx=x[edgeidx[1][i]]-x[edgeidx[0][i]];
        double dy=y[edgeidx[1][i]]-y[edgeidx[0][i]];
        barx_[0][nBarx_*i  ]=hatx_[0][i]-dx/c;
        barx_[1][nBarx_*i  ]=hatx_[1][i]-dy/c;
        barx_[0][nBarx_*i+1]=hatx_[0][i]+dx/c;
        barx_[1][nBarx_*i+1]=hatx_[1][i]+dy/c;
        edge_len_[i]=sqrt(dx*dx+dy*dy);
    }

    // 计算周长
    for (unsigned long i=0; i<nElement; ++i){
        perimeter_[i] = 0.0;
        for(int j=0; j<nSide_; j++){  // 三角形三条边
        perimeter_[i] += edge_len_[tri_edge_conn[j][i]];
        }
    }
}


// 计算外法向量
void LinearDGSolver_2D::computeOuterNormal(){
    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    unsigned long ** tri=trimesh_->triangle();
    unsigned long ** tri_edge_conn = trimesh_->tri_edge_connection();

    for(unsigned long ik=0; ik<nElement; ik++){
        for(int ie=0; ie<nSide_; ie++){  
            // 该边的两个点索引
            unsigned long p1 = tri[ie][ik];
            unsigned long p2 = tri[(ie < 2)? ie + 1 : 0][ik];
            // 边编号
            unsigned long idx_edge = tri_edge_conn[ie][ik];
            // idx_edge边的两个点
            double vx1 = x[p1];
            double vy1 = y[p1];
            double vx2 = x[p2];
            double vy2 = y[p2];
            outerNormal_[0][ik*nSide_+ie] = (vy2 - vy1) / edge_len_[idx_edge];
            outerNormal_[1][ik*nSide_+ie] = (vx1 - vx2) / edge_len_[idx_edge];
        }
    }
}


// 计算基函数在相关节点上的值
void LinearDGSolver_2D::computeBasisOnNodes(){
    unsigned long ** tri_edge_conn = trimesh_->tri_edge_connection();

    // 计算基函数相关值 
    // 在不同单元顶点基函数值是相同的
    phi_vec_[0][0] =-1.0; phi_vec_[0][1] = 1.0; phi_vec_[0][2] = 1.0;
    phi_vec_[1][0] = 1.0; phi_vec_[1][1] =-1.0; phi_vec_[1][2] = 1.0;
    phi_vec_[2][0] = 1.0; phi_vec_[2][1] = 1.0; phi_vec_[2][2] =-1.0;
    for(unsigned long ik=0; ik<nElement; ik++){ //按单元循环
        for(int ie=0; ie<nSide_; ie++){ //按边循环
            double phival[nPhi_], phival_x[nPhi_], phival_y[nPhi_];
            unsigned long edgeIndex = tri_edge_conn[ie][ik];
            phi(ik, hatx_[0][edgeIndex], hatx_[1][edgeIndex], phival, phival_x, phival_y);
            for(int j=0; j<nPhi_; j++){
                phi_hatx_x_[j][ik*nHatx_+ie] = phival_x[j];  // nHatx_ = nSide_ 
                phi_hatx_y_[j][ik*nHatx_+ie] = phival_y[j];
            }
            for(int m=0; m<nBarx_; m++){  // 每条边上的两个高斯点
                phi(ik, barx_[0][nBarx_*edgeIndex+m], barx_[1][nBarx_*edgeIndex+m], phival);
                for(int j=0; j<nPhi_; j++){
                    phi_barx_[j][nBarx_*(ik*nSide_+ie)+m] = phival[j]; 
                }
            }
        }
    }
}

// 计算边界信息
void LinearDGSolver_2D::computeBoundaryInfo(){

    unsigned long ** edgeidx=trimesh_->edge_info();
    // 记录边界个数和边界的边索引
    nBoundry = 0;
    unsigned long tempboundry[nEdge];
    for(unsigned long i=0; i < nEdge; i++){
        if(edgeidx[3][i] == nElement+1){
            tempboundry[nBoundry] = i;
            nBoundry++;
        }                
    }
    boundary_ = new unsigned long[nBoundry];
    for(unsigned long i=0; i < nBoundry; i++){
        boundary_[i] = tempboundry[i]; 
    }  
    
}



// 计算初始状态
void LinearDGSolver_2D::computeInitialData() {
    unsigned long** conn=trimesh_->tri_edge_connection();
    double *u_tmp[nVars_];
    for (int i=0; i<nVars_; i++){
        u_tmp[i]=new double[nEdge];
        u0_[i] = new double[dof_];
    }

    // 初值条件
    computeInitialCondition(nEdge, hatx_[0], hatx_[1], u_tmp);

    // 赋值给u_
    for (unsigned long i = 0; i < nElement ; i++){
        for (int j = 0; j < nPhi_ ; j++) {
            for (int k=0; k<nVars_; k++){
                u_[k][i*nPhi_+j]=u_tmp[k][conn[j][i]];
                u0_[k][i*nPhi_+j]=u_tmp[k][conn[j][i]];
            }
        }
    }
    for (int i=0; i<nVars_; i++){
        delete []u_tmp[i];
    }
}


void LinearDGSolver_2D::computeFluxOnElement(unsigned long ik, double** u, double f[2][12]){
    double fu[dim_][4]; // 用于储存该边的f(u)通量
    double ui[nVars_];  // 储存该边的物理量
    for(int j=0; j<nPhi_; j++){  // 三个基函数
        for(int n=0; n<nVars_; n++){
            ui[n] = u[n][ik*nPhi_+j];  // ui为hatx(中点)上的数值解 由于采用中点插值 所以数值解即等于系数
        }
        flux(ui, fu);  // 计算该边通量
        // 将该边通量储存到f中 f即储存了一个单元三条边的通量
        for(int k=0; k<dim_; k++){
            for(int n=0; n<nVars_; n++){
                f[k][j*nVars_+n] = fu[k][n];
            }
        }
    }
}


void LinearDGSolver_2D::computeNumericalSolutionOnGaussPoints(unsigned long ik, unsigned long ie, double** u, double uh[2][4]){
    for (int i=0; i<nVars_; i++){      // 遍历4个物理量
        for (int m=0; m<nBarx_; m++){
            uh[m][i]=0;
            for (int j=0; j<nPhi_; j++){  // 遍历基函数
                uh[m][i] += u[i][ik*nPhi_+j] * phi_barx_[j][nBarx_*(ik*nSide_+ie)+m];
            }
        }
    }
}



// 计算单元内部流通量, 返回2x4矩阵
void LinearDGSolver_2D::computeNumericalFluxOnInteriorEdge(unsigned long ik, unsigned long ie, double** u, double f[2][4]){
    
    unsigned long** edgeInfo = trimesh_->edge_info();
    unsigned long** tri_edge_conn = trimesh_->tri_edge_connection();
    unsigned long** edge_order_in_tri = trimesh_->edge_order_in_tri();

    // 边编号
    unsigned long idx_edge = tri_edge_conn[ie][ik];    

    // 外法向量
    double normal[dim_];
    normal[0] = outerNormal_[0][ik*nSide_+ie];
    normal[1] = outerNormal_[1][ik*nSide_+ie];

    // 该边相邻的两个三角形
    unsigned long k1 = edgeInfo[2][idx_edge];
    unsigned long k2 = edgeInfo[3][idx_edge];

    // idx_edge对应的另一个三角形ik_ 以及 idx_edge在ik_中的order ie_
    unsigned long ik_, ie_; 
    (k1 == ik) ? (ik_ = k2, ie_ = edge_order_in_tri[1][idx_edge]) : 
    (k2 == ik) ? (ik_ = k1, ie_ = edge_order_in_tri[0][idx_edge]) :
    (ik_ = nElement+1, ie_ = 3);
    if (ik_ == nElement+1) {  // 找不到 (边界)
        std::cerr << "Error: can't find another triangle" << std::endl;
        throw -1;
    }

    // 计算u(x, y) 在高斯点上的值
    double u1[nBarx_][4], u2[nBarx_][4];
    computeNumericalSolutionOnGaussPoints(ik, ie, u, u1);
    computeNumericalSolutionOnGaussPoints(ik_, ie_, u, u2);

    // 计算流通量
    for (int m=0; m<nBarx_; m++){  // 对于两个高斯点
        double ff[nVars_];
        numerical_flux(u1[m], u2[m], normal, ff);
        for(int n=0; n<nVars_; n++){ // 对于4个物理量
            f[m][n] = ff[n]; 
        }
    }
}
    



// 计算整个单元边界上的数值流通量, 返回6x4矩阵
void LinearDGSolver_2D::computeNumericalFluxOnElement(unsigned long ik, double** u, double f[6][4]){

    unsigned long** edgeInfo = trimesh_->edge_info();
    unsigned long **tri_edge_conn_ = trimesh_->tri_edge_connection();
    
    double flux[dim_][4];
    for(int i=0; i<nSide_; i++){  // 三条边
        unsigned long ie = tri_edge_conn_[i][ik];  // 该边索引
        // 判断该边是否是边界
        if (edgeInfo[3][ie] == nElement+1){  // 是边界
            computeNumericalFluxOnBoundary(ik, i, u, flux);

        } 
        else{
            computeNumericalFluxOnInteriorEdge(ik, i, u, flux);
        }
        
        // 将该边通量储存到f中 f即储存了一个单元三条边(6个点)的数值通量
        for(int k=0; k<nBarx_; k++){
            for(int j=0; j<nVars_; j++){
                f[nBarx_*i+k][j] = flux[k][j];
            }
        }
    }
}


// 计算空间离散, 返回 4x dof_ 的矩阵f
void LinearDGSolver_2D::computeSpaceDiscretization(double** u, double**f) {

    unsigned long** conn=trimesh_->tri_edge_connection();
    for(int n=0; n<nVars_; n++){
        for (unsigned long j=0; j<dof_; j++){
            f[n][j]=0.0;
        }
    }
    double fu[dim_][12], flux[nBarx_*nSide_][4];
    for (unsigned long ik=0; ik<nElement; ik++){  // 单元数
        computeFluxOnElement(ik, u, fu);  // 计算k单元的f(u)通量
        computeNumericalFluxOnElement(ik, u, flux);  // 计算k单元的数值通量
        double omege = area_[ik]/3;
        for (int j=0; j<nPhi_; j++){  // 基函数循环
            // 常微分方程的第一项
            for (int i=0; i<nHatx_; i++){ // 积分点循环
                for(int n=0; n<nVars_; n++){   // 对于4个物理量
                    f[n][ik*nPhi_+j] += fu[0][i*nVars_+n] * phi_hatx_x_[j][ik*nHatx_+i] + fu[1][i*nVars_+n] * phi_hatx_y_[j][ik*nHatx_+i]; 
                }   
            }
            // 常微分方程的第二项
            for (int ie=0; ie<nSide_; ie++){   // 边循环
                double scale2 =  edge_len_[conn[ie][ik]]/(nBarx_*omege);
                for (int m=0; m<nBarx_; m++){  // 边上高斯点循环
                    for(int n=0; n<nVars_; n++){  // 对于4个物理量循环
                        f[n][ik*nPhi_+j] -= scale2 * flux[nBarx_*ie+m][n] * phi_barx_[j][nBarx_*(ik*nSide_+ie)+m];
                    }
                }
            }
        }
    }
}

double LinearDGSolver_2D::computeTimeStep(double** u){
    
    // 条件数
    double CFL = 0.3;
    double dm[nElement]={0}; 

    for(unsigned long i=0; i<nElement; i++){
        double u_average[nVars_];
        for(int j=0; j<nVars_; j++){
            u_average[j] = (u[j][i*nPhi_] + u[j][i*nPhi_+1] + u[j][i*nPhi_+2]) / nPhi_;
        }
        double rho=u_average[0], E=u_average[3], v[dim_]={u_average[1]/u_average[0], u_average[2]/u_average[0]};
        double cs =SOUND_SPEED(PRESSURE(E, rho, v), rho);
        dm[i] = (l2_norm(v, dim_) + cs) * (perimeter_[i] / area_[i]);
    }

    int n = sizeof(dm) / sizeof(dm[0]); // 数组大小 (nElement)
    double max_dm = *std::max_element(dm, dm + n);  // 寻找最大值 

    return CFL/max_dm;
    
}


void LinearDGSolver_2D::computeTimeDiscretization(double total_time, double* save_time_points, int save_time_num, double dt){
    
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
void LinearDGSolver_2D::NumericalSolutionOnVertex(double** u, int idx, double* var){
    unsigned long **conn=trimesh_->triangle();
    
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
void LinearDGSolver_2D::NumericalSolutionOnVertex(double** u, double** var){

    unsigned long **conn=trimesh_->triangle();

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
void LinearDGSolver_2D::SynchronizationUpdate(double** u){
    NumericalSolutionOnVertex(u, 0, rho);
}

    

void LinearDGSolver_2D::outputTecPlotDataFile(const std::string &fname, double time_point){

    double *x=trimesh_->x_coord();
    double *y=trimesh_->y_coord();
    unsigned long **conn=trimesh_->triangle();

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
    output << "TITLE=\"Triangular mesh\"" << std::endl;
    output << "VARIABLES= \"x\" \"y\" \"rho\" \"vx\" \"vy\" \"E\" \"p\"" << std::endl;
    output << "ZONE T=\"none\"," << "N=" << nVertex << "," << "E=" << nElement << "," << "ET=TRIANGLE, F=FEPOINT" << std::endl;
    
    for (unsigned long i = 0; i < nVertex; i++) {
        double rho = uh[0][i];
        double vx = uh[1][i] / rho;
        double vy = uh[2][i] / rho;
        double v[dim_] = {vx, vy};
        double E = uh[3][i];
        double p = PRESSURE(E, rho, v);

        output << std::setprecision(15) << std::setw(20) << x[i] << "\t";
        output << std::setprecision(15) << std::setw(20) << y[i] << "\t";
        output << std::setprecision(15) << std::setw(20) << rho << "\t";
        output << std::setprecision(15) << std::setw(20) << vx << "\t";
        output << std::setprecision(15) << std::setw(20) << vy << "\t";
        output << std::setprecision(15) << std::setw(20) << E << "\t";
        output << std::setprecision(15) << std::setw(20) << p << std::endl;
    }

    for (unsigned long i = 0; i < nElement; i++) {
        output << std::setw(10) << conn[0][i] + 1 << "\t";
        output << std::setw(10) << conn[1][i] + 1 << "\t";
        output << std::setw(10) << conn[2][i] + 1 << std::endl;
    }
    output.close();

    std::cout << "successfully written to " << fname << std::endl;

    for (int i = 0; i < nVars_; i++){
        delete[] uh[i];
        delete[] uj[i];
    }
}

void LinearDGSolver_2D::outputVTKDataFile(const std::string &fname, double time_point) {
    double *x = trimesh_->x_coord();
    double *y = trimesh_->y_coord();
    unsigned long **conn = trimesh_->triangle();

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
        double v[dim_] = {vx, vy};
        double E = uh[3][i];
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

