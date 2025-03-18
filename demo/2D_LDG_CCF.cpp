#include <iostream>
#include "LinearDGSolver2DCavityCylinderFlow.h"
#include <filesystem> 

using namespace std;

int main() {
    
    std::string read_path = std::string(SOURCE_DIR) + "/data/2D/CavityCylinderFlow/";
    std::string write_path = std::string(SOURCE_DIR) + "/result/2D/CavityCylinderFlow/";
    std::string fileName = "rectangle_with_center_circle0.1.off";

    TriangleMesh *mesh=new TriangleMesh;
    mesh->read_off(read_path + fileName);
    mesh->collect_edges();

    LinearDGSolver_2D* dle2Solver = new LinearDGSolver_2D_CavityCylinderFlow(mesh);

    double total_time = 3.5;

    int save_time_num = 3;
    double save_time[save_time_num] = {0, 1, 3.5};

    dle2Solver->computeTimeDiscretization(total_time, save_time, save_time_num);

    string filename0 = "0s_result.dat";
    string filename1 = "1s_result.dat";           
    string filename2 = "3.5s_result.dat";                                                                                  
    dle2Solver->outputTecPlotDataFile(write_path + filename0, 0);
    dle2Solver->outputTecPlotDataFile(write_path + filename1, 1);
    dle2Solver->outputTecPlotDataFile(write_path + filename2, 3.5);

    delete mesh;
    delete dle2Solver;


    return 0;
}