#include "LinearDGSolver2DCavityCylinderFlow.h"

// 计算初始状态
void LinearDGSolver_2D_CavityCylinderFlow::computeInitialCondition(unsigned long n, double *x, double *y, double **rtval) {
    for (unsigned long i=0; i<n; i++){
        if(x[i]<=-1){
            rtval[0][i] = 1.4764;
            rtval[1][i] = 0.4877;
            rtval[2][i] = 0.0;
            rtval[3][i] = 1.7372;
        }
        if(x[i]>-1){
            rtval[0][i] = 1.0;
            rtval[1][i] = 0.0;
            rtval[2][i] = 0.0;
            rtval[3][i] = 1.0;
        }
    }
}