#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_LineSearchStep.hpp>

#include "vector.hpp"

#include <fstream>
#include <cereal/archives/binary.hpp>

void serialise_step(const ROL::Step<double> &v) {
  std::ofstream os("rol_step.cereal", std::ios::binary);
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(v);
}

void load_step(ROL::Step<double> &v) {
  std::ifstream is("rol_step.cereal", std::ios::binary);
  cereal::BinaryInputArchive iarchive(is);

  iarchive(v);
}

void init_step(py::module& m) {
  m.def("serialise_step", &serialise_step, "Serialise a ROL Step using Cereal");
  m.def("load_step", &load_step, "Load a ROL Step from Cereal archive");

  py::class_<ROL::Step<double>, std::shared_ptr<ROL::Step<double>>>(m, "Step");

  py::class_<ROL::LineSearchStep<double>, std::shared_ptr<ROL::LineSearchStep<double>>, ROL::Step<double>>(m, "LineSearchStep")
    .def(py::init([](ROL::ParameterList &params, std::shared_ptr<ROL::Secant<double>> &secant) {
      return new ROL::LineSearchStep<double>(params, nullptr, secant);
    }));
}
