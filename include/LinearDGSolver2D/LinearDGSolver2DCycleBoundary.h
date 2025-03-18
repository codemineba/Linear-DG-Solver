#ifndef _LINEAR_DG_SOLVER_2D_CYCLE_BOUNDARY_H_
#define _LINEAR_DG_SOLVER_2D_CYCLE_BOUNDARY_H_


#include "LinearDGSolver2D.h"
#define PI acos(-1)

class LinearDGSolver_2D_CycleBoundary: public LinearDGSolver_2D {
    
private:
    unsigned long *periodic_boundary_;  // 用于储存周期边界信息

public: 

    LinearDGSolver_2D_CycleBoundary(TriangleMesh *mesh)
        : LinearDGSolver_2D(mesh) {
            // 计算周期边界信息
            computePeriodicBoundaryInfo();
        }

    ~LinearDGSolver_2D_CycleBoundary() override{
        delete[] periodic_boundary_;
    }
    
    // 计算周期边界信息
    void computePeriodicBoundaryInfo();

    // 计算初始状态
    void computeInitialCondition(unsigned long n, double *x, double *y, double **rtval) override ;

    // 计算边界边上的数值通量
    void computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[2][4]) override ;

    // 计算rho的误差的L2范数
    double computeL2ErrorOfRho();

    // 输出日志信息
    std::string LogMessage();

};



#endif // _LINEAR_DG_SOLVER_2D_CYCLE_BOUNDARY_H_