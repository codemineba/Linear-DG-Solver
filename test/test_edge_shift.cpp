#include <iostream>
#include <cmath>

using namespace std;

// 用于判断两个 double 是否近似相等
bool is_almost_equal(double a, double b, double scale = 1e-9) {
    return fabs(a - b) < scale;
}

// 判断 edge1 是否可以通过平移得到 edge2
bool is_edge_shift_from(double edge1[2][2], double edge2[2][2]) {
    // 检查 edge1 和 edge2 是否为水平边
    bool is_edge1_horizontal = is_almost_equal(edge1[0][1], edge1[1][1]);
    bool is_edge2_horizontal = is_almost_equal(edge2[0][1], edge2[1][1]);

    // 检查 edge1 和 edge2 是否为竖直边
    bool is_edge1_vertical = is_almost_equal(edge1[0][0], edge1[1][0]);
    bool is_edge2_vertical = is_almost_equal(edge2[0][0], edge2[1][0]);

    // 检查是否完全重合
    bool is_completely_overlapping = 
        is_almost_equal(edge1[0][0], edge2[0][0]) &&
        is_almost_equal(edge1[0][1], edge2[0][1]) &&
        is_almost_equal(edge1[1][0], edge2[1][0]) &&
        is_almost_equal(edge1[1][1], edge2[1][1]);

    // 如果完全重合，则返回 false
    if (is_completely_overlapping) {
        return false;
    }

    // 如果 edge1 和 edge2 同为水平边
    if (is_edge1_horizontal && is_edge2_horizontal) {
        // 检查两个点的 x 坐标是否完全相等
        return is_almost_equal(edge1[0][0], edge2[0][0]) &&
               is_almost_equal(edge1[1][0], edge2[1][0]);
    }

    // 如果 edge1 和 edge2 同为竖直边
    if (is_edge1_vertical && is_edge2_vertical) {
        // 检查两个点的 y 坐标是否完全相等
        return is_almost_equal(edge1[0][1], edge2[0][1]) &&
               is_almost_equal(edge1[1][1], edge2[1][1]);
    }

    // 不满足水平或竖直的条件
    return false;
}


int main() {
    // 测试边集合
    double edge1[2][2] = {{0, 1}, {2, 1}};  // 水平边
    double edge2[2][2] = {{0, 1}, {2, 1}};  // 完全重合
    double edge3[2][2] = {{0, 2}, {2, 2}};  // 水平边 y 平移
    double edge4[2][2] = {{1, 0}, {1, 2}};  // 竖直边
    double edge5[2][2] = {{1, 1}, {1, 3}};  // 竖直边 y 平移
    double edge6[2][2] = {{1, 1}, {3, 1}};  // 水平边，但 x 坐标不同
    double edge7[2][2] = {{1, 0}, {1, 2}};  // 与 edge4 相同
    double edge8[2][2] = {{2, 1}, {4, 1}};  // 平移得到新的水平边
    double edge9[2][2] = {{0, 10.0/1214}, {2, 10.0/1214}};  // 竖直边，但 x 不同

    // 测试多个边组合
    cout << "Testing edges:" << endl;
    cout << "edge1 vs edge2: " << (is_edge_shift_from(edge1, edge2) ? "True" : "False") << endl;  // False
    cout << "edge1 vs edge3: " << (is_edge_shift_from(edge1, edge3) ? "True" : "False") << endl;  // True
    cout << "edge1 vs edge4: " << (is_edge_shift_from(edge1, edge4) ? "True" : "False") << endl;  // False
    cout << "edge4 vs edge5: " << (is_edge_shift_from(edge4, edge5) ? "True" : "False") << endl;  // False
    cout << "edge1 vs edge6: " << (is_edge_shift_from(edge1, edge6) ? "True" : "False") << endl;  // False
    cout << "edge4 vs edge7: " << (is_edge_shift_from(edge4, edge7) ? "True" : "False") << endl;  // False
    cout << "edge1 vs edge8: " << (is_edge_shift_from(edge1, edge8) ? "True" : "False") << endl;  // False
    cout << "edge1 vs edge9: " << (is_edge_shift_from(edge1, edge9) ? "True" : "False") << endl;  // True

    return 0;
}