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
	 py::arg("obj"), py::arg("x"), py::arg("g") = (ROL::Vector<double>*)nullptr)
    .def("addBoundConstraint", &ROL::Problem<double>::addBoundConstraint)
    .def("addInequalityConstraint",
	 py::overload_cast<
	 std::string,
	 const ROL::Ptr<ROL::Constraint<double>>&,
	 const ROL::Ptr<ROL::Vector<double>>&,
	 const ROL::Ptr<ROL::BoundConstraint<double>>&,
	 const ROL::Ptr<ROL::Vector<double>>&,
	 bool>(&ROL::Problem<double>::addConstraint),
	 py::arg("name"), py::arg("icon"), py::arg("imul"), py::arg("ibnd"),
	 py::arg("ires") = (ROL::Vector<double>*)nullptr,
	 py::arg("reset") = false)
    .def("finalize", [](ROL::Problem<double> &instance) {
      instance.finalize(false, true);
    });
}
