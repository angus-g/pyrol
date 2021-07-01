#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <map>

#include <ROL_Problem.hpp>

void init_problem(py::module &m) {
  py::class_<ROL::Problem<double>,
	     std::shared_ptr<ROL::Problem<double>>>(m, "Problem")
    .def(py::init<std::shared_ptr<ROL::Objective<double>>&,
	 ROL::Ptr<ROL::Vector<double>>&,
	 ROL::Ptr<ROL::Vector<double>>&>(),
	 py::arg("obj"), py::arg("x"), py::arg("g") = (ROL::Vector<double>*)nullptr);
}
