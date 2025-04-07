#include <pybind11/pybind11.h>
#include <pybind11/functional.h>  
#include <pybind11/numpy.h>
#include <functional>
#include <memory>
#include "LinearDGSolver2DCycleBoundary.h"
#include "LinearDGSolver2DCavityCylinderFlow.h"

namespace py = pybind11;


// 返回值是一维数组
py::array_t<double> getRenderingData1(LinearDGSolver_2D_CycleBoundary& self) {
    // 获取 density_ 数组的指针
    double* arrayPtr = self.getRenderingData();

    // 获取数组的大小
    size_t size = self.getNRenderingData();

    // 创建 NumPy 数组，使用 buffer_info 来管理内存
    auto array = py::array(py::buffer_info(
        arrayPtr,             // 指向数据的指针
        sizeof(double),       // Size of one scalar
        py::format_descriptor<double>::format(), // Python struct-style format descriptor
        1,                    // Number of dimensions
        { size },             // Shape of the array
        { sizeof(double) }    // Strides (in bytes) for each axis
    ));

    // 返回 NumPy 数组
    return array;
}


// 解决Python GIL线程锁释放的问题
void computeTimeDiscretization1(LinearDGSolver_2D_CycleBoundary& self, double total_time) {
    // Release GIL to allow other Python threads to run
    py::gil_scoped_release release;
    
    self.computeTimeDiscretization(total_time);
    
    // Re-acquire GIL before returning to Python
    py::gil_scoped_acquire acquire;
}



// 返回值是一维数组
py::array_t<double> getRenderingData2(LinearDGSolver_2D_CavityCylinderFlow& self) {
    // 获取 density_ 数组的指针
    double* arrayPtr = self.getRenderingData();

    // 获取数组的大小
    size_t size = self.getNRenderingData();

    // 创建 NumPy 数组，使用 buffer_info 来管理内存
    auto array = py::array(py::buffer_info(
        arrayPtr,             // 指向数据的指针
        sizeof(double),       // Size of one scalar
        py::format_descriptor<double>::format(), // Python struct-style format descriptor
        1,                    // Number of dimensions
        { size },             // Shape of the array
        { sizeof(double) }    // Strides (in bytes) for each axis
    ));

    // 返回 NumPy 数组
    return array;
}


// 解决Python GIL线程锁释放的问题
void computeTimeDiscretization2(LinearDGSolver_2D_CavityCylinderFlow& self, double total_time) {
    // Release GIL to allow other Python threads to run
    py::gil_scoped_release release;
    
    self.computeTimeDiscretization(total_time);
    
    // Re-acquire GIL before returning to Python
    py::gil_scoped_acquire acquire;
}




PYBIND11_MODULE(PY_DG_2D, m)
{

    py::class_<TriangleMesh>(m, "TriangleMesh")
        .def(py::init<>())
        .def("read_off", &TriangleMesh::read_off)
        .def("collect_edges", &TriangleMesh::collect_edges);

    py::class_<LinearDGSolver_2D_CycleBoundary>(m, "LinearDGSolver_2D_CycleBoundary")
        .def(py::init<TriangleMesh*>())
        .def("getRenderingData", &getRenderingData1) 
        .def("computeTimeDiscretization", &computeTimeDiscretization1);

    py::class_<LinearDGSolver_2D_CavityCylinderFlow>(m, "LinearDGSolver_2D_CavityCylinderFlow")
        .def(py::init<TriangleMesh*>())
        .def("getRenderingData", &getRenderingData2) 
        .def("computeTimeDiscretization", &computeTimeDiscretization2);

}
