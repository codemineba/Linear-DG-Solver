#include "mesh_structure/TriangleMesh.h"


// 标准C++不提供去除首尾空格的库函数，但利用string也可以实现此功能
std::string TriangleMesh::trim(std::string s) {
    if (s.empty()) { return s; }
    s.erase(0, s.find_first_not_of(' '));
    s.erase(s.find_last_not_of(' ') + 1);
    // 注意linux和Windows文件格式的区别：
    // 当linux上的代码读取Windows文件格式时， 读取结果的每行都会多一个\r
    // 因而此处需要删除\r
    size_t pos=s.find_last_of('\r');    // 当linux上的代码读取Windows文件格式时， 读取结果的每行都会多一个\r
    if (pos != std::string::npos){      // 因而可能需要删除这个额外的\r
      s.erase(pos);
    }   
    return s;
}

bool TriangleMesh::splitOff(const std::string &str, const std::string &pattern, int npattern, std::string *res)
{
    // 先去除首尾空格
    std::string text = trim(str);
    if(str.empty()) return false;
    //在字符串末尾加入分隔符，是为了截取最后一段
    std::string strs = text + pattern;
    size_t pos = strs.find(pattern);  // size_t提高程序的可移植性，可以理解为无符号int
    int n = 0;
    // string::npos为一个静态成员常量，表示size_t的最大值，通常用来表示"直到字符串结尾"
    while(pos != std::string::npos)   // 如果还有索引位置
    {
        std::string temp = strs.substr(0, pos);
        res[n++] = temp;
        //去掉已分割的字符串,在剩下的字符串中进行分割
        strs = strs.substr(pos + 1, strs.size());
        pos = strs.find(pattern);
        if (n >= npattern)  break;     // 增加容错性
    }
    return true;
}


void TriangleMesh::read_off(const std::string &path) {
    std::fstream file(path, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Wrong file name or path!" << std::endl;
        throw -1;
    }

    std::string line;
    getline(file, line);
    if (trim(line) != "OFF") {
        std::cerr << "Not a valid OFF header!" << std::endl;
        throw -1;
    }

    getline(file, line);
    std::string first_line[4];
    if (splitOff(line, " ", 2, first_line)){
        nVertex   = stoi(first_line[0]);
        nTriangle = stoi(first_line[1]);
    }
    x_ = new double [nVertex];
    y_ = new double [nVertex];
    z_ = new double [nVertex];
    
    tri_[0] = new unsigned long[nTriangle];
    tri_[1] = new unsigned long[nTriangle];
    tri_[2] = new unsigned long[nTriangle];

    std::string data[4];
    for (unsigned long i=0; i<nVertex; ++i){
        getline(file, line);
        if (splitOff(line," ", 3, data)){
            x_[i] = stod(data[0]);
            y_[i] = stod(data[1]);
            z_[i] = stod(data[2]);
        }   
    }

    for (unsigned long i=0; i<nTriangle; ++i){
        getline(file, line);
        if(splitOff(line," ", 4, data)){
            tri_[0][i] = stoi(data[1]);
            tri_[1][i] = stoi(data[2]);
            tri_[2][i] = stoi(data[3]);

            // 规定三角形索引必须按照逆时针排序
            // 计算该三角形的法向量
            double normal[3];
            check_point_order(i, normal);
            // 如果三角形是顺时针的，则交换两个点，调整为逆时针 
            if (normal[2] < 0) { 
                std::swap(tri_[1][i], tri_[2][i]);
            }
        }
    }
}



// 计算两点的叉积，用于判断三角形法向量的方向
void TriangleMesh::check_point_order(unsigned long ik, double* normal) {
    
    double x1 = x_[tri_[0][ik]], y1 = y_[tri_[0][ik]], z1 = z_[tri_[0][ik]], 
           x2 = x_[tri_[1][ik]], y2 = y_[tri_[1][ik]], z2 = z_[tri_[1][ik]],
           x3 = x_[tri_[2][ik]], y3 = y_[tri_[2][ik]], z3 = z_[tri_[2][ik]];

    double u[3] = {x2-x1, y2-y1, z2-z1};
    double v[3] = {x3-x1, y3-y1, z3-z1};
    
    // 计算normal叉积
    normal[0] = u[1] * v[2] - u[2] * v[1];
    normal[1] = u[2] * v[0] - u[0] * v[2];
    normal[2] = u[0] * v[1] - u[1] * v[0];
}


// 寻找与点相邻的单元
void TriangleMesh::find_vertex_trangle_connection(std::vector<unsigned long> *conn) {
    for (unsigned long i =0; i<nTriangle; ++i){
        for (int j=0; j<3; j++){
            conn[tri_[j][i]].push_back(i);
        }
    }
}



void TriangleMesh::collect_edges(){
    int max_nedge_each_vertex = 100;  
    std::vector<std::tuple<unsigned long, unsigned long, unsigned long, unsigned long>> edges;
    int *localedge[max_nedge_each_vertex];
    nEdge = 0;
    for (int i = 0; i < max_nedge_each_vertex; i++){
        localedge[i] = new int[nVertex];
    }
    int *nlocaledge = new int[nVertex];
    for (unsigned long j = 0; j < nVertex; j++){
        nlocaledge[j] = 0;
    }

    tri_edge_conn_[0] = new unsigned long[nTriangle];
    tri_edge_conn_[1] = new unsigned long[nTriangle];
    tri_edge_conn_[2] = new unsigned long[nTriangle];

    for (unsigned long i = 0; i < nTriangle; ++i){   // 按单元循环
        for (int j = 0; j < 3; j++){
            unsigned long start, end;  // end > start
            start = tri_[j][i] > tri_[(j + 1) % 3][i] ? tri_[(j + 1) % 3][i] : tri_[j][i];
            end = tri_[j][i] < tri_[(j + 1) % 3][i] ? tri_[(j + 1) % 3][i] : tri_[j][i];
            int k = 0;
            for (k = 0; k < nlocaledge[start]; k++){    // 判断边是否已经存在
                std::tuple<unsigned long, unsigned long, unsigned long, unsigned long> &t = edges[localedge[k][start]];
                if (std::get<1>(t) == end){       // 边已经存在
                    std::get<3>(t) = i;
                    tri_edge_conn_[j][i] = localedge[k][start];
                    break;
                }
            }
            if (k == nlocaledge[start]){          // 边不存在
                std::tuple<unsigned long, unsigned long, unsigned long, unsigned long> t(start, end, i, nTriangle+1);
                edges.push_back(t); 
                localedge[nlocaledge[start]][start] = nEdge;
                tri_edge_conn_[j][i] = localedge[k][start];
                ++nEdge;
                ++nlocaledge[start];    
                assert(nlocaledge[start] < max_nedge_each_vertex);
            }
        }
    }

    for (int i = 0; i < max_nedge_each_vertex; i++){
        delete[] localedge[i];
    }
    delete[] nlocaledge;

    for (int i = 0; i < 4; i++){
        edges_[i] = new unsigned long[nEdge];
    }
    for (unsigned long j = 0; j < nEdge; ++j){
        std::tuple<unsigned long, unsigned long, unsigned long, long> t = edges[j];
        edges_[0][j] = std::get<0>(t);
        edges_[1][j] = std::get<1>(t);
        edges_[2][j] = std::get<2>(t);
        edges_[3][j] = std::get<3>(t);

    }
    for (int i = 0; i < 2; i++){
        edge_order_in_tri_[i] = new unsigned long[nEdge];
        for(unsigned long j=0; j<nEdge; j++) {edge_order_in_tri_[i][j] = 3;}
    }
    for (unsigned long j = 0; j < nEdge; ++j){
        unsigned long ik = edges_[2][j], ik_ = edges_[3][j];
        for (unsigned long i = 0; i < 3; ++i){  // order of ik and ik_
            if (tri_edge_conn_[i][ik] == j) { edge_order_in_tri_[0][j] = i;}
            if (tri_edge_conn_[i][ik_] == j) {edge_order_in_tri_[1][j] = i;}
        } 
    }
}


void TriangleMesh::outputTecPlotDataFile(const char *fname){
	std::ofstream output(fname);
    output<<"TITLE=\"Triangular mesh\""<<std::endl;
    output<<"VARIABLES= \"x\" \"y\" \"z\" \"var\""<<std::endl;
    output<<"ZONE T=\"none\","<<"N="<<nVertex<<","<<"E="<<nTriangle<<","<<"ET=TRIANGLE, F=FEPOINT"<<std::endl;
    for (unsigned long i=0; i<nVertex; i++){
        output<<std::setprecision(15)<<std::setw(20)<<x_[i]<<"\t";
        output<<std::setprecision(15)<<std::setw(20)<<y_[i]<<"\t";
        output<<std::setprecision(15)<<std::setw(20)<<z_[i]<<"\t";
        output<<std::setprecision(15)<<std::setw(20)<<0<<std::endl;
    }
    
    for (unsigned long i=0; i<nTriangle; i++){
        output<<std::setw(10)<<tri_[0][i]<<"\t";
        output<<std::setw(10)<<tri_[1][i]<<"\t";
        output<<std::setw(10)<<tri_[2][i]<<std::endl;
    }
	output.close();
}