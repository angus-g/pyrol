#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_TypeB_LinMoreAlgorithm.hpp>
#include <ROL_TypeB_QuasiNewtonAlgorithm.hpp>
#include <ROL_Bounds.hpp>

void init_algorithm(py::module &m) {
  py::class_<ROL::TypeB::LinMoreAlgorithm<double>>(m, "LinMoreAlgorithm")
    .def(py::init<ROL::ParameterList&>())
    .def(py::init<ROL::ParameterList&, const std::shared_ptr<ROL::Secant<double>>&>())
    .def("setStatusTest", &ROL::TypeB::LinMoreAlgorithm<double>::setStatusTest)
    .def("run",
	 [](ROL::TypeB::LinMoreAlgorithm<double> &instance, ROL::Vector<double> &x,
	    ROL::Objective<double> &obj, ROL::Bounds<double> &bnd) {
	   instance.run(x, obj, bnd, std::cout);
	 });

  py::class_<ROL::TypeB::QuasiNewtonAlgorithm<double>>(m, "QuasiNewtonAlgorithm")
    .def(py::init<ROL::ParameterList&>())
    .def(py::init<ROL::ParameterList&, const std::shared_ptr<ROL::Secant<double>>&>())
    .def("setStatusTest", &ROL::TypeB::QuasiNewtonAlgorithm<double>::setStatusTest)
    .def("run",
	 [](ROL::TypeB::QuasiNewtonAlgorithm<double> &instance, ROL::Vector<double> &x,
	    ROL::Objective<double> &obj, ROL::Bounds<double> &bnd) {
	   instance.run(x, obj, bnd,std::cout);
	 });
}
