#include <iostream>
#include <filesystem> 
#include "LinearDGSolver3DCycleBoundary.h"
#include "mesh_structure/TetrahedronMesh.h"

using namespace std;

int main() {

    std::string read_path = std::string(SOURCE_DIR) + "/data/3D/CycleBoundary/";
    std::string write_path = std::string(SOURCE_DIR) + "/result/3D/CycleBoundary/";
    std::string fileName = "tetrahedronMesh0.2.off";
    
    TetrahedronMesh *mesh=new TetrahedronMesh;
    mesh->read_off(read_path + fileName);    
    mesh->collect_faces();

    LinearDGSolver_3D* dle3Solver = new LinearDGSolver_3D_CycleBoundary(mesh);

    double total_time = 2;

    int save_time_num = 1;
    double save_time[save_time_num] = {0};
    // double dt = 0.001;

    dle3Solver->computeTimeDiscretization(total_time, save_time, save_time_num);

    // string filename0 = "0s_result.dat";
    // string filename1 = "1s_result.dat";                                                                                        
    // string filename2 = "2s_result.dat";
    // dle3Solver->outputTecPlotDataFile(write_path + filename0, 0);
    // dle3Solver->outputTecPlotDataFile(write_path + filename1, 1);
    // dle3Solver->outputTecPlotDataFile(write_path + filename2, 2);


    delete mesh;
    delete dle3Solver;

    return 0;
}