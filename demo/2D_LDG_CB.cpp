#include <iostream>
#include <filesystem>
#include "LinearDGSolver2DCycleBoundary.h" 

using namespace std;

int main() {

    std::string read_path = std::string(SOURCE_DIR) + "/data/2D/CycleBoundary/";
    std::string write_path = std::string(SOURCE_DIR) + "/result/2D/CycleBoundary/";
    std::string fileName = "triangleMesh0.0125.off";
    
    TriangleMesh *mesh=new TriangleMesh;
    mesh->read_off(read_path + fileName);
    mesh->collect_edges();

    LinearDGSolver_2D* dle2Solver = new LinearDGSolver_2D_CycleBoundary(mesh);

    double total_time = 2;

    int save_time_num = 3;
    double save_time[save_time_num] = {0, 1, 2};
    // double dt = 0.001;

    dle2Solver->computeTimeDiscretization(total_time, save_time, save_time_num);

    string filename0 = "0s_result.dat";
    string filename1 = "1s_result.dat";                                                                                        
    string filename2 = "2s_result.dat";
    dle2Solver->outputTecPlotDataFile(write_path + filename0, 0);
    dle2Solver->outputTecPlotDataFile(write_path + filename1, 1);
    dle2Solver->outputTecPlotDataFile(write_path + filename2, 2);


    delete mesh;
    delete dle2Solver;

    return 0;
}