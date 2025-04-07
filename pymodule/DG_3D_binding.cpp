#include <pybind11/pybind11.h>
#include <pybind11/functional.h>  
#include <pybind11/numpy.h>
#include <functional>
#include <memory>
#include "LinearDGSolver3DCycleBoundary.h"

namespace py = pybind11;


// 返回值是一维数组
py::array_t<double> getRenderingData(LinearDGSolver_3D_CycleBoundary& self) {
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
void computeTimeDiscretization(LinearDGSolver_3D_CycleBoundary& self, double total_time) {
    // Release GIL to allow other Python threads to run
    py::gil_scoped_release release;
    
    self.computeTimeDiscretization(total_time);
    
    // Re-acquire GIL before returning to Python
    py::gil_scoped_acquire acquire;
}




PYBIND11_MODULE(PY_DG_3D, m)
{

    py::class_<TetrahedronMesh>(m, "TetrahedronMesh")
        .def(py::init<>())
        .def("read_off", &TetrahedronMesh::read_off)
        .def("collect_faces", &TetrahedronMesh::collect_faces);

    py::class_<LinearDGSolver_3D_CycleBoundary>(m, "LinearDGSolver_3D_CycleBoundary")
        .def(py::init<TetrahedronMesh*>())
        .def("getRenderingData", &getRenderingData) 
        .def("computeTimeDiscretization", &computeTimeDiscretization);
}
