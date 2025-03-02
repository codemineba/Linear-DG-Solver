#ifndef MESH_STRUCTURE_TEXT_TETRAHEDRONMESH_H
#define MESH_STRUCTURE_TEXT_TETRAHEDRONMESH_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <tuple>
#include <cassert>
#include <algorithm>

using namespace std;

class TetrahedronMesh {

private:
    unsigned long nVertex;       // 顶点数
    unsigned long nEdge;         // 边数
    unsigned long nFace;         // 面数
    unsigned long nTetrahedron;  // 四面体数
    unsigned long nBoundary;     // 边界数

    double *x_;
    double *y_;
    double *z_;
    unsigned long *boundary_;

    unsigned long *tet_[4];
    unsigned long *tet_face_conn_[4];
    unsigned long *faces_[5];  // index of 3 vertices and 2 neighboring tetrahedron for each edge
    unsigned long *face_order_in_tet_[2];    // order of edge in 2 neighboring tetrahedron


    // 标准C++不提供去除首尾空格的库函数，但利用string也可以实现此功能
    string trim(string s);

    bool splitOff(const string &str, const string &pattern, int n, string *res);

    void check_point_order(unsigned long ik, double* var);  // 用于判断点的排序

public:
    TetrahedronMesh(){
        nVertex=0;
        nEdge=0;
        nFace=0;
        nTetrahedron=0;
        nBoundary=0;

        x_=nullptr;
        y_=nullptr;
        z_=nullptr;

        tet_[0]=tet_[1]=tet_[2]=tet_[3]=nullptr;
        tet_face_conn_[0]=tet_face_conn_[1]=tet_face_conn_[2]=tet_face_conn_[3]=nullptr;
        faces_[0]=faces_[1]=faces_[2]=faces_[3]=faces_[4]=nullptr;
        face_order_in_tet_[0]=face_order_in_tet_[1]=nullptr;
    }

    ~TetrahedronMesh(){
        if (nTetrahedron!=0){
            delete []x_;
            delete []y_;
            delete []z_;
            delete []boundary_;
            delete []tet_[0];
            delete []tet_[1];
            delete []tet_[2];
            delete []tet_[3];
        }
        if (nFace!=0){
            for (int i=0; i<5; i++){
                delete []faces_[i];
            }
            for (int i=0; i<4; i++){
                delete []tet_face_conn_[i];
            }

            delete []face_order_in_tet_[0];
            delete []face_order_in_tet_[1];
        }
        nVertex=0;
        nEdge=0;
        nFace=0;
        nTetrahedron=0;
    }

    void read_off(const string &path) ;

    unsigned long getNVertex(){  return nVertex; }
    unsigned long getNFace(){  return nFace; }
    unsigned long getNTetrahedron() { return nTetrahedron; }
    unsigned long getNBoundary() { return nBoundary; }
    double *x_coord(){  return x_; }
    double *y_coord(){  return y_; }
    double *z_coord(){  return z_; }
    unsigned long *boundary(){ return boundary_; }

    unsigned long **tetrahedron()  { return tet_; }
    unsigned long **face_info() { return faces_; }
    unsigned long **tet_face_connection() { return tet_face_conn_; }
    unsigned long **face_order_in_tet() { return face_order_in_tet_; }

    void find_vertex_tetrahedron_connection(std::vector<unsigned long> *conn);

    void collect_faces();

    void outputTecPlotDataFile(const char *fname);
};    // TetrahedronMesh

#endif //MESH_STRUCTURE_TEXT_TETRAHEDRONMESH_H
