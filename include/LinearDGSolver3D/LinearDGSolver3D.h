#ifndef _LINEAR_DG_SOLVER_2D_H_
#define _LINEAR_DG_SOLVER_2D_H_


#include <cmath>
#include <cstring> 
#include <algorithm>
#include "mesh_structure/TetrahedronMesh.h"
#include "SciCalUtils/GeometryUtils.h"

// 定义宏
#define GAMMA 1.4  // 绝热系数
#define CFL 0.3    // 条件数
#define PRESSURE(E, rho, v) ((GAMMA - 1) * ((E) - 0.5 * (rho) * (v[0]*v[0] + v[1]*v[1]+v[2]*v[2])))  // 压强
#define SOUND_SPEED(p, rho) (sqrt(GAMMA * p / rho))  // 声速


class LinearDGSolver_3D {

protected:

    const int dim_;                // 维数
    const int nVars_;              // 物理量(守恒量)个数
    const int nHatx_;              // 单元内求积点个数(每个单元)  
    const int nBarx_;              // 单元边界上求积点个数(每条边)
    const int nPhi_;               // 基函数个数
    const int nSide_;              // 单元的边(面)数
    const int nCorner_;            // 单元的角(点)数

    unsigned long dof_;         // 单元个数乘4, 每个变量的整体自由度
    double *u_[5];              // 解向量u
    double *u0_[5];             // 记录初始时刻的值
    double *hatx_[3];           // 储存标准单元内的求积节点
    double *barx_[3];           // 储存标准单元边界上的求积节点
    double *phi_hatx_[4];       // hatx的基函数值
    double *phi_hatx_x_[4];     // hatx的基函数x方向导数值
    double *phi_hatx_y_[4];     // hatx的基函数y方向导数值
    double *phi_hatx_z_[4];     // hatx的基函数z方向导数值
    double *phi_barx_[4];       // barx的基函数值
    double *phi_vec_[4];        // 顶点基函数值
    double *outerNormal_[3];    // 储存所有单元的外法向量
    double *face_area_;         // 面的面积，用于计算外法向量
    double *sufacearea_;        // 单元表面积
    double *volume_;            // 单元体积 
    double *rho;


    unsigned long nElement;     // 单元数
    unsigned long nVertex;      // 顶点数
    unsigned long nFace;        // 面数
    unsigned long nBoundary;     // 边界数
    
    double timeInterval_;          // 时间步长    
    double timePoint_;             // 实时时刻 
    double *times_point_save_;     // 需要保存的时间点(默认只保存开始时间和终止时间的计算结果)
    int num_times_point_save_;     // 需要保存的时间点的个数

    TetrahedronMesh *tetmesh_;


public:

    // 构造函数
    LinearDGSolver_3D(TetrahedronMesh *mesh) 
        : dim_(3), nVars_(5), nHatx_(4), nBarx_(3), nPhi_(4), nSide_(4) ,nCorner_(4) {

        timeInterval_ = 0.0;
        timePoint_ = 0.0;
        
        tetmesh_= mesh;
        nElement = mesh->getNTetrahedron();
        nVertex = mesh->getNVertex();
        nFace = mesh->getNFace();
        nBoundary = mesh->getNBoundary();
        dof_=nElement*nPhi_;

        // 数组初始化
        face_area_= new double[nFace];         // 面积
        volume_ = new double[nElement];        // 体积
        sufacearea_ = new double[nElement];    // 表面积
        for (int i=0; i<dim_; i++){ 
            hatx_[i] = new double[nElement*nHatx_];   // 标准单元内求积节点
            barx_[i] = new double[nFace*nBarx_];      // 标准单元边界上求积节点   
            outerNormal_[i] = new double[dof_];  // 外法向量
        }
        for(int i=0; i<nPhi_; i++){
            phi_hatx_[i] = new double[dof_];  // hatx的基函数x方向导数值
            phi_hatx_x_[i] = new double[dof_];  // hatx的基函数x方向导数值
            phi_hatx_y_[i] = new double[dof_];  // hatx的基函数y方向导数值
            phi_hatx_z_[i] = new double[dof_];  // hatx的基函数z方向导数值
            phi_barx_[i] = new double[dof_*nBarx_];  // barx的基函数值
            phi_vec_[i] = new double[nCorner_];        // 顶点基函数值
        }
        rho = new double[nVertex];  // 为绘制实时渲染记录变量
        for(unsigned long i=0; i<nVertex; i++){
            rho[i] = 0.0;
        }
        computeElementProperties();    // 计算单元属性 (包括: 内部求积点,边界求积点,面积,表面积,体积)
        computeOuterNormal();          // 计算外法向量
        computeBasisOnNodes();         // 计算基函数在相关节点上的值

    }
    
    // 析构函数
    virtual ~LinearDGSolver_3D() {
        delete[] volume_;
        delete[] face_area_;
        delete[] sufacearea_;
        for (int i=0; i<3; i++){
            delete[] hatx_[i];
            delete[] barx_[i];
        }
        for (int i=0; i<3; i++){ 
            delete[] outerNormal_[i];
        }
        for (int i=0; i<4; i++){
            delete[] phi_hatx_x_[i];
            delete[] phi_hatx_y_[i];
            delete[] phi_hatx_z_[i];
            delete[] phi_barx_[i];
            delete[] phi_vec_[i];
        }
        for (int i=0; i<5; i++){
            delete[] u_[i];
            delete[] u0_[i];
        }
        delete[] rho;
    }


    double* getVar() {return rho;}
    void setVar(double* var) { rho = var;}

    // 基函数及其导数
    void phi(unsigned long ik, double x, double y, double z, double *var, double *var_x=nullptr, double *var_y=nullptr, double *var_z=nullptr);

    // 通量
    void flux(double* u, double f[3][5]);

    // 数值通量
    void numerical_flux(double* u1, double *u2, double *normal, double* f);

    // 计算单元属性  (包括: 内部求积点,边界求积点,面积,表面积,体积)
    void computeElementProperties();

    // 计算基函数在相关节点上的值
    void computeBasisOnNodes();
    
    // 计算外法向量
    void computeOuterNormal();

    // 计算初始状态
    void computeInitialData();

    // 计算初值条件
    virtual void computeInitialCondition(unsigned long n, double *x, double *y, double *z, double **rtval) = 0;

    // 计算边界边上的数值通量
    virtual void computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[3][5]) = 0;

    // 计算内部边上数值通量
    void computeNumericalFluxOnInteriorEdge(unsigned long ik, unsigned long ie, double** u, double f[3][5]);

    // 计算单元边界上的数值通量
    void computeNumericalFluxOnElement(unsigned long ik, double** u, double f[12][5]);

    // 计算单元边界上的通量
    void computeFluxOnElement(unsigned long ik, double** u, double f[3][20]);

    // 计算空间离散
    void computeSpaceDiscretization(double** u, double** f);

    // 计算时间离散
    void computeTimeDiscretization(double total_time, double* save_time_points=nullptr, int save_time_num=2, double dt=0.0);

    // 计算某个时刻的单元顶点的四个数值解
    void NumericalSolutionOnVertex(double** u, int idx, double* var);

    // 计算某个时刻的单元顶点的全部数值解
    void NumericalSolutionOnVertex(double** u, double** var);

    // 计算边上高斯点的数值解
    void computeNumericalSolutionOnGaussPoints(unsigned long ik, unsigned long ie, double** u, double uh[3][5]);

    // 计算时间步长
    double computeTimeStep(double** u);

    // 日志信息
    virtual std::string LogMessage() =0 ;

    // 将结果保存至tecplot数据格式的文件
    void outputTecPlotDataFile(const std::string &fname, double time_point);

    // 将结果保存至VTK数据格式的文件
    void outputVTKDataFile(const std::string &fname, double time_point);

    // 同步更新解值
    void SynchronizationUpdate(double** u);
 
};

#endif // _LINEAR_DG_SOLVER_2D_H_
