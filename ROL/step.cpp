#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_LineSearchStep.hpp>
#include <ROL_TrustRegionStep.hpp>

void init_step(py::module& m) {
  py::class_<ROL::Step<double>, std::shared_ptr<ROL::Step<double>>>(m, "Step");

  py::class_<ROL::LineSearchStep<double>, std::shared_ptr<ROL::LineSearchStep<double>>, ROL::Step<double>>(m, "LineSearchStep")
    .def(py::init([](ROL::ParameterList &params, std::shared_ptr<ROL::Secant<double>> &secant) {
      return new ROL::LineSearchStep<double>(params, nullptr, secant);
    }));

  py::class_<ROL::TrustRegionStep<double>, std::shared_ptr<ROL::TrustRegionStep<double>>, ROL::Step<double>>(m, "TrustRegionStep")
    .def(py::init<std::shared_ptr<ROL::Secant<double>>&, ROL::ParameterList&>());
}
