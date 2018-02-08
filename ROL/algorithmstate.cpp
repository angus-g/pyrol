#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Types.hpp>

void init_algorithmstate(py::module& m) {
    py::class_<ROL::AlgorithmState<double>,
               std::shared_ptr<ROL::AlgorithmState<double>>>(m,
                                                             "AlgorithmState")
        .def_readwrite("gnorm", &ROL::AlgorithmState<double>::gnorm)
        .def_readwrite("cnorm", &ROL::AlgorithmState<double>::cnorm)
        .def_readwrite("snorm", &ROL::AlgorithmState<double>::snorm)
        .def_readwrite("value", &ROL::AlgorithmState<double>::value);
}
