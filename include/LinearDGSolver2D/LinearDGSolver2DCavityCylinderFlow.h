#ifndef _LINEAR_DG_SOLVER_2D_CAVITY_CYLINDER_FLOW_H_
#define _LINEAR_DG_SOLVER_2D_CAVITY_CYLINDER_FLOW_H_

#include "LinearDGSolver2DCavityObstacleFlow.h"

class LinearDGSolver_2D_CavityCylinderFlow: public LinearDGSolver_2D_CavityObstacleFlow {

public: 

    LinearDGSolver_2D_CavityCylinderFlow(TriangleMesh *mesh)
        : LinearDGSolver_2D_CavityObstacleFlow(mesh){
            top_boundary_value_=8.0;       // 上边界
            bottom_boundary_value_=-8.0;   // 下边界
            left_boundary_value_=-8.0;     // 左边界
            right_boundary_value_=8.0;     // 右边界
            // 判断方腔内的边界类型
            computeBoundaryType();
        }

    ~LinearDGSolver_2D_CavityCylinderFlow() override  = default;

    // 计算初始状态
    void computeInitialCondition(unsigned long n, double *x, double *y, double **rtval) override ;
};


#endif  // _LINEAR_DG_SOLVER_2D_CAVITY_CYLINDER_FLOW_H_
