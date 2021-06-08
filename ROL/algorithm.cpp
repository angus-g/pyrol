#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "vector.hpp"

#include <ROL_Algorithm.hpp>
#include <ROL_Bounds.hpp>

#include <fstream>
#include <cereal/archives/binary.hpp>

void serialise_algorithm(const ROL::Algorithm<double> &v) {
  std::ofstream os("rol_algorithm.cereal", std::ios::binary);
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(v);
}

void load_algorithm(ROL::Algorithm<double> &v) {
  std::ifstream is("rol_algorithm.cereal", std::ios::binary);
  cereal::BinaryInputArchive iarchive(is);

  iarchive(v);
}

void init_algorithm(py::module& m) {
  m.def("serialise_algorithm", &serialise_algorithm, "Serialise a ROL Algorithm using Cereal");
  m.def("load_algorithm", &load_algorithm, "Load a ROL Algorithm from Cereal archive");

  //
  // ROL::Algorithm<double>
  //
  py::class_<ROL::Algorithm<double>>(m, "Algorithm")
    .def(py::init<const std::string&, ROL::ParameterList&>())
    .def(py::init<const std::shared_ptr<ROL::Step<double>>&, const std::shared_ptr<ROL::StatusTest<double>>&>())
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Objective<double>& obj) {
	   instance.run(x, obj, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Objective<double>& obj, ROL::Bounds<double>& bnd) {
	   instance.run(x, obj, bnd, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Vector<double>& l, ROL::Objective<double>& obj,
	    ROL::Constraint<double>& con) {
	   instance.run(x, l, obj, con, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
	    ROL::Vector<double>& l, ROL::Objective<double>& obj,
	    ROL::Constraint<double>& con, ROL::Bounds<double>& bnd) {
	   instance.run(x, l, obj, con, bnd, true, std::cout);
	 })
    .def("run",
	 [](ROL::Algorithm<double>& instance,
	    ROL::OptimizationProblem<double>& opt) {
	   instance.run(opt, true, std::cout);
	 })
    .def("get_state", [](ROL::Algorithm<double>& instance) {
      return instance.getState();
    });
}
