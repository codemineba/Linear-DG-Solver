#include "mesh_structure/TetrahedronMesh.h"

// 标准C++不提供去除首尾空格的库函数，但利用string也可以实现此功能
string TetrahedronMesh::trim(string s) {
    if (s.empty()) { return s; }
    s.erase(0, s.find_first_not_of(' '));
    s.erase(s.find_last_not_of(' ') + 1);
    // 注意linux和Windows文件格式的区别：
    // 当linux上的代码读取Windows文件格式时， 读取结果的每行都会多一个\r
    // 因而此处需要删除\r
    size_t pos=s.find_last_of('\r');    // 当linux上的代码读取Windows文件格式时， 读取结果的每行都会多一个\r
    if (pos != std::string::npos){         // 因而可能需要删除这个额外的\r
        s.erase(pos);
    }
    return s;
}

bool TetrahedronMesh::splitOff(const string &str, const string &pattern, int npattern, string *res)
{
    // 先去除首尾空格
    string text = trim(str);
    if(str.empty()) return false;
    //在字符串末尾加入分隔符，是为了截取最后一段
    string strs = text + pattern;
    size_t pos = strs.find(pattern);  // size_t提高程序的可移植性，可以理解为无符号int
    int n = 0;
    // string::npos为一个静态成员常量，表示size_t的最大值，通常用来表示"直到字符串结尾"
    while(pos != string::npos)   // 如果还有索引位置
    {
        string temp = strs.substr(0, pos);
        res[n++] = temp;
        //去掉已分割的字符串,在剩下的字符串中进行分割
        strs = strs.substr(pos + 1, strs.size());
        pos = strs.find(pattern);
        if (n >= npattern)  break;     // 增加容错性
    }
    return true;
}


void TetrahedronMesh::read_off(const string &path) {
    fstream file(path, ios::in);
    if (!file.is_open()) {
        cout << "Wrong file name or path!" << endl;
        throw -1;
    }

    string line;
    getline(file, line);
    if (trim(line) != "OFF") {
        cout << "Not a valid OFF header!" << endl;
        throw -1;
    }

    getline(file, line);
    string first_line[4];
    if (splitOff(line, " ", 2, first_line)){
        nVertex   = stoi(first_line[0]);
        nTetrahedron = stoi(first_line[1]);
    }
    x_ = new double [nVertex];
    y_ = new double [nVertex];
    z_ = new double [nVertex];

    tet_[0] = new unsigned long[nTetrahedron];
    tet_[1] = new unsigned long[nTetrahedron];
    tet_[2] = new unsigned long[nTetrahedron];
    tet_[3] = new unsigned long[nTetrahedron];

    string data[5];
    for (unsigned long i=0; i<nVertex; ++i){
        getline(file, line);
        if (splitOff(line," ", 3, data)){
            x_[i] = stod(data[0]);
            y_[i] = stod(data[1]);
            z_[i] = stod(data[2]);
        }
    }

    for (unsigned long i=0; i<nTetrahedron; ++i){
        getline(file, line);
        if(splitOff(line," ", 5, data)){
            tet_[0][i] = stoi(data[1]);
            tet_[1][i] = stoi(data[2]);
            tet_[2][i] = stoi(data[3]);
            tet_[3][i] = stoi(data[4]);
        }
    }
}

// 用于判断四面体点的排序
// void TetrahedronMesh::check_point_order(unsigned long ik, double* var) {
    
//     double x0 = x_[tet_[0][ik]], y0 = y_[tet_[0][ik]], z0 = z_[tet_[0][ik]],
//            x1 = x_[tet_[1][ik]], y1 = y_[tet_[1][ik]], z1 = z_[tet_[1][ik]], 
//            x2 = x_[tet_[2][ik]], y2 = y_[tet_[2][ik]], z2 = z_[tet_[2][ik]],
//            x3 = x_[tet_[3][ik]], y3 = y_[tet_[3][ik]], z3 = z_[tet_[3][ik]];
           
//     double v0[3]={x0, y0, z0};
//     double v1[3]={x1, y1, z1};
//     double v2[3]={x2, y2, z2};
//     double v3[3]={x3, y3, z3};

// }

// 寻找与点相邻的单元
void TetrahedronMesh::find_vertex_tetrahedron_connection(std::vector<unsigned long> *conn){
    for (unsigned long i =0; i< nTetrahedron; ++i){
        for (int j=0; j<4; j++){
            conn[tet_[j][i]].push_back(i);
        }
    }
}

void TetrahedronMesh::collect_faces(){
    int max_nface_each_vertex = 300;
    std::vector<std::tuple<unsigned long, unsigned long,unsigned long,unsigned long,unsigned long>> faces;
    int** localface = new int* [max_nface_each_vertex];
    nFace=0;
    for (int i=0; i<max_nface_each_vertex; i++){
        localface[i]=new int[nVertex];
    }
    int *nlocalface=new int[nVertex];
    for (unsigned long j=0; j<nVertex; j++){
        nlocalface[j]=0;
    }

    tet_face_conn_[0]=new unsigned long[nTetrahedron];
    tet_face_conn_[1]=new unsigned long[nTetrahedron];
    tet_face_conn_[2]=new unsigned long[nTetrahedron];
    tet_face_conn_[3]=new unsigned long[nTetrahedron];

    int iv[4][3]={{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {3, 2, 1}};
    for (unsigned long i =0; i<nTetrahedron; ++i){   // 按单元循环
        for (int j=0; j<4; j++){  //每个单元的4个面
            std::vector<unsigned long> giv= {tet_[iv[j][0]][i], tet_[iv[j][1]][i], tet_[iv[j][2]][i]};
            std::sort(giv.begin(), giv.end());
            int k;
            for (k=0; k<nlocalface[giv[0]]; k++){    // 判断面是否已经存在
                std::tuple<unsigned long, unsigned long,unsigned long,unsigned long,unsigned long> &f=faces[localface[k][giv[0]]];
                if (std::get<1>(f)==giv[1] && std::get<2>(f)==giv[2]){
                    std::get<4>(f)=i;     // 面已经存在
                    tet_face_conn_[j][i]=localface[k][giv[0]];
                    break;
                }
            }
            if (k==nlocalface[giv[0]]){          // 面不存在
                std::tuple<unsigned long, unsigned long,unsigned long,unsigned long,unsigned long> t(giv[0], giv[1], giv[2], i, nTetrahedron+1);

                faces.push_back(t);
                localface[nlocalface[giv[0]]][giv[0]]=nFace;
                tet_face_conn_[j][i]=localface[k][giv[0]];
                ++nFace;
                ++nlocalface[giv[0]];
                assert(nlocalface[giv[0]]< max_nface_each_vertex);
            }
        }
    }


    for (int i=0; i<max_nface_each_vertex; i++){
        delete []localface[i];
    }
    delete []nlocalface;

    for (int i=0; i<5; i++){
        faces_[i]=new unsigned long [nFace];
    }
    for (unsigned long j=0; j<nFace; ++j){
        std::tuple<unsigned long, unsigned long, unsigned long, unsigned long, unsigned long> t=faces[j];
        faces_[0][j]=std::get<0>(t);
        faces_[1][j]=std::get<1>(t);
        faces_[2][j]=std::get<2>(t);
        faces_[3][j]=std::get<3>(t);
        faces_[4][j]=std::get<4>(t);
    }
    for (int i = 0; i < 2; i++){
        face_order_in_tet_[i] = new unsigned long[nFace];
        for(unsigned long j=0; j<nFace; j++) {face_order_in_tet_[i][j] = 4;}
    }
    for (unsigned long j = 0; j < nFace; ++j){
        unsigned long ik = faces_[3][j], ik_ = faces_[4][j];
        for (unsigned long i = 0; i < 4; ++i){  // order of ik and ik_
            if (tet_face_conn_[i][ik] == j) { face_order_in_tet_[0][j] = i;}
            if (tet_face_conn_[i][ik_] == j) {face_order_in_tet_[1][j] = i;}
        } 
    }

}


void TetrahedronMesh::outputTecPlotDataFile(const char *fname){
    ofstream output(fname);
    output<<"TITLE=\"Quadrilateral mesh\""<<endl;
    output<<"VARIABLES= \"x\" \"y\" \"z\" \"var\""<<endl;
    output<<"ZONE T=\"none\","<<"N="<<nVertex<<","<<"E="<<nTetrahedron<<","<<"ET=TETRAHEDRON, F=FEPOINT"<<endl;
    for (unsigned long i=0; i<nVertex; i++){
        output<<setprecision(15)<<setw(20)<<x_[i]<<"\t";
        output<<setprecision(15)<<setw(20)<<y_[i]<<"\t";
        output<<setprecision(15)<<setw(20)<<z_[i]<<"\t";
        output<<setprecision(15)<<setw(20)<<0<<endl;
    }

    for (unsigned long i=0; i<nTetrahedron; i++){
        output<<setw(10)<<tet_[0][i]+1<<"\t";
        output<<setw(10)<<tet_[1][i]+1<<"\t";
        output<<setw(10)<<tet_[2][i]+1<<"\t";
        output<<setw(10)<<tet_[3][i]+1<<endl;
    }
    output.close();
}