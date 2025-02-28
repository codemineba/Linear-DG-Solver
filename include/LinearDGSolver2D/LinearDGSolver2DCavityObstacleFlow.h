#ifndef _LINEAR_DG_SOLVER_2D_CAVITY_OBSTACLE_FLOW_H_
#define _LINEAR_DG_SOLVER_2D_CAVITY_OBSTACLE_FLOW_H_

#include "LinearDGSolver2D.h"

class LinearDGSolver_2D_CavityObstacleFlow: public LinearDGSolver_2D {
    
protected:

    double top_boundary_value_;     // 上边界
    double bottom_boundary_value_;  // 下边界
    double left_boundary_value_;    // 左边界
    double right_boundary_value_;   // 右边界

    unsigned long *boundary_type_;  // 用于储存边界类型 0: 入流边界 1: 出流边界 2: 反射边界
    
public: 

    LinearDGSolver_2D_CavityObstacleFlow(TriangleMesh *mesh)
        : LinearDGSolver_2D(mesh){}

    ~LinearDGSolver_2D_CavityObstacleFlow() override{
        delete boundary_type_;
    }

    // 计算边界类型
    void computeBoundaryType();

    // 计算边界边上的数值通量
    void computeNumericalFluxOnBoundary(unsigned long ik, unsigned long ie, double** u, double f[2][4]) override ;

    // 计算Rho的2范数
    double computeL2NormOfRho();

    // 输出日志信息
    std::string LogMessage();

};



#endif // _LINEAR_DG_SOLVER_2D_CAVITY_OBSTACLE_FLOW_H_