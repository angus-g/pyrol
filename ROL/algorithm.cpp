#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include "vector.hpp"

#include <ROL_Algorithm.hpp>
#include <ROL_TypeB_LinMoreAlgorithm.hpp>
#include <ROL_TypeB_QuasiNewtonAlgorithm.hpp>
#include <ROL_Bounds.hpp>

#include <fstream>
#include <filesystem>
#define FMT_HEADER_ONLY
#include <fmt/core.h>

CEREAL_REGISTER_TYPE(ROL::TypeB::LinMoreAlgorithm<double>);

void serialise_algorithm(const std::shared_ptr<ROL::TypeB::Algorithm<double>> &v, int rank, std::string checkpoint_dir) {
  std::filesystem::path checkpoint_path(checkpoint_dir);
  std::ofstream os(checkpoint_path.append(fmt::format("rol_typeb_algorithm_{}.cereal", rank)), std::ios::binary);
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(v);
}

std::shared_ptr<ROL::TypeB::Algorithm<double>> load_algorithm(int rank, std::string checkpoint_dir) {
  std::filesystem::path checkpoint_path(checkpoint_dir);
  std::ifstream is(checkpoint_path.append(fmt::format("rol_typeb_algorithm_{}.cereal", rank)), std::ios::binary);
  cereal::BinaryInputArchive iarchive(is);

  std::shared_ptr<ROL::TypeB::Algorithm<double>> v;
  iarchive(v);
  return v;
}

void init_algorithm(py::module& m) {
  m.def("serialise_algorithm", &serialise_algorithm, "Serialise a ROL Algorithm using Cereal",
	py::arg("algorithm"), py::arg("rank") = 0, py::arg("checkpoint_dir") = "");
  m.def("load_algorithm", &load_algorithm, "Load a ROL Algorithm from Cereal archive",
	py::arg("rank") = 0, py::arg("checkpoint_dir") = "");

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

  py::class_<ROL::TypeB::Algorithm<double>, std::shared_ptr<ROL::TypeB::Algorithm<double>>>(m, "TypeBAlgorithm");

  // ROL 2.0 TypeB algorithms
  py::class_<ROL::TypeB::LinMoreAlgorithm<double>, ROL::TypeB::Algorithm<double>, std::shared_ptr<ROL::TypeB::LinMoreAlgorithm<double>>>(m, "LinMoreAlgorithm")
    .def(py::init<ROL::ParameterList&>())
    .def(py::init<ROL::ParameterList&, const std::shared_ptr<ROL::Secant<double>>&>())
    .def("setStatusTest", &ROL::TypeB::LinMoreAlgorithm<double>::setStatusTest)
    .def("run",
	 [](ROL::TypeB::LinMoreAlgorithm<double> &instance, ROL::Vector<double> &x,
	    ROL::Objective<double> &obj, ROL::Bounds<double> &bnd) {
	   instance.run(x, obj, bnd, std::cout);
	 });
}
