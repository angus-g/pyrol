#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>

#include "vector.hpp"

#include <ROL_Algorithm.hpp>
#include <ROL_TypeB_LinMoreAlgorithm.hpp>
#include <ROL_TypeB_QuasiNewtonAlgorithm.hpp>
#include <ROL_Bounds.hpp>

#include <fstream>
#define FMT_HEADER_ONLY
#include <fmt/core.h>

void serialise_algorithm(const ROL::TypeB::Algorithm<double> &v, int rank) {
  std::ofstream os(fmt::format("rol_typeb_algorithm_{}.cereal", rank), std::ios::binary);
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(v);
}

void load_algorithm(ROL::TypeB::Algorithm<double> &v, int rank) {
  std::ifstream is(fmt::format("rol_typeb_algorithm_{}.cereal", rank), std::ios::binary);
  cereal::BinaryInputArchive iarchive(is);

  iarchive(v);
}

void init_algorithm(py::module& m) {
  m.def("serialise_algorithm", &serialise_algorithm, "Serialise a ROL Algorithm using Cereal",
	py::arg("algorithm"), py::arg("rank") = 0);
  m.def("load_algorithm", &load_algorithm, "Load a ROL Algorithm from Cereal archive",
	py::arg("algorithm"), py::arg("rank") = 0);

  //
  // ROL::Algorithm<double>
  //
  py::class_<ROL::Algorithm<double>>(m, "Algorithm")
    .def(py::init<const std::shared_ptr<ROL::Step<double>>&, const std::shared_ptr<ROL::StatusTest<double>>&>())
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Objective<double>& obj) {
	   // type U interface
	   instance.run(x, obj, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Objective<double>& obj, ROL::Bounds<double>& bnd) {
	   // type B interface
	   instance.run(x, obj, bnd, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Vector<double>& l, ROL::Objective<double>& obj,
	    ROL::Constraint<double>& con) {
	   // type E interface
	   instance.run(x, l, obj, con, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Vector<double>& l, ROL::Objective<double>& obj,
	    ROL::Constraint<double>& con, ROL::Bounds<double>& bnd) {
	 // type EB interface
	 instance.run(x, l, obj, con, bnd, true, std::cout);
	 });

  py::class_<ROL::TypeB::Algorithm<double>>(m, "TypeBAlgorithm");

  // ROL 2.0 TypeB algorithms
  py::class_<ROL::TypeB::LinMoreAlgorithm<double>, ROL::TypeB::Algorithm<double>>(m, "LinMoreAlgorithm")
    .def(py::init<ROL::ParameterList&>())
    .def(py::init<ROL::ParameterList&, const std::shared_ptr<ROL::Secant<double>>&>())
    .def("setStatusTest", &ROL::TypeB::LinMoreAlgorithm<double>::setStatusTest)
    .def("run",
	 [](ROL::TypeB::LinMoreAlgorithm<double> &instance, ROL::Vector<double> &x,
	    ROL::Objective<double> &obj, ROL::Bounds<double> &bnd) {
	   instance.run(x, obj, bnd, std::cout);
	 });

  py::class_<ROL::TypeB::QuasiNewtonAlgorithm<double>, ROL::TypeB::Algorithm<double>>(m, "QuasiNewtonAlgorithm")
    .def(py::init<ROL::ParameterList&>())
    .def(py::init<ROL::ParameterList&, const std::shared_ptr<ROL::Secant<double>>&>())
    .def("setStatusTest", &ROL::TypeB::QuasiNewtonAlgorithm<double>::setStatusTest)
    .def("run",
	 [](ROL::TypeB::QuasiNewtonAlgorithm<double> &instance, ROL::Vector<double> &x,
	    ROL::Objective<double> &obj, ROL::Bounds<double> &bnd) {
	   instance.run(x, obj, bnd,std::cout);
	 });
}
