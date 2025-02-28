#ifndef _TRIANGLE_MESH_H_
#define _TRIANGLE_MESH_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <tuple>
#include <assert.h>

class TriangleMesh {  // 三角形的点和边和面的索引都是从0开始

private:
    unsigned long nVertex;
    unsigned long nEdge;
    unsigned long nTriangle;

    double *x_;
    double *y_;
    double *z_;

    unsigned long *tri_[3];
    unsigned long *tri_edge_conn_[3];
    unsigned long *edges_[4];  // index of 2 vertices and 2 neighboring triangle for each edge 
    unsigned long *edge_order_in_tri_[2];    // order of edge in 2 neighboring triangle

    // 标准C++不提供去除首尾空格的库函数，但利用string也可以实现此功能
    std::string trim(std::string s);

    bool splitOff(const std::string &str, const std::string &pattern, int n, std::string *res);

    void check_point_order(unsigned long ik, double* normal);  // 计算ik单元的外法向量 用于判断点的排序

public:
    TriangleMesh(){
        nVertex=0;
        nTriangle=0;

        x_=nullptr;
        y_=nullptr;
        z_=nullptr;

        tri_[0]=tri_[1]=tri_[2]=nullptr;
        tri_edge_conn_[0]=tri_edge_conn_[1]=tri_edge_conn_[2]=nullptr;
        edges_[0]=edges_[1]=edges_[2]=edges_[3]=nullptr;
        edge_order_in_tri_[0]=edge_order_in_tri_[1]=nullptr;
    }
    
    ~TriangleMesh(){
        if (nTriangle!=0){
            delete []x_;
            delete []y_;
            delete []z_;
            delete []tri_[0];
            delete []tri_[1];
            delete []tri_[2];
        }
        if (nEdge!=0){
            delete []edges_[0];
            delete []edges_[1];
            delete []edges_[2];
            delete []edges_[3];

            delete []tri_edge_conn_[0];
            delete []tri_edge_conn_[1];
            delete []tri_edge_conn_[2];
            
            delete []edge_order_in_tri_[0];
            delete []edge_order_in_tri_[1];
        }
    }

    void read_off(const std::string &path) ;

    unsigned long getNVertex(){  return nVertex; }
    unsigned long getNEdge()  {  return nEdge;   }
    unsigned long getNTriangle() { return nTriangle; }
    double *x_coord(){  return x_; }
    double *y_coord(){  return y_; }
    double *z_coord(){  return z_; }

    unsigned long **triangle()  { return tri_; }
    unsigned long **edge_info() { return edges_;}
    unsigned long **tri_edge_connection() {  return tri_edge_conn_; }
    unsigned long **edge_order_in_tri() { return edge_order_in_tri_; }

    void collect_edges();

    void find_vertex_trangle_connection(std::vector<unsigned long> *conn);
    
    void outputTecPlotDataFile(const char *fname);
};    // TriangleMesh

#endif
