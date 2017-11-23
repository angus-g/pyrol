#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

#include <ROL_OptimizationProblem.hpp>
#include <ROL_OptimizationSolver.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include "Teuchos_ParameterList.hpp"

namespace py = pybind11;
void init_optimizationsolver(py::module& m) {
  // ROL::OptimizationSolver<double>
  //

  py::class_<ROL::OptimizationSolver<double>,
             std::shared_ptr<ROL::OptimizationSolver<double>>>(
      m, "OptimizationSolver")
      .def("__init__",
           [](ROL::OptimizationSolver<double>& instance,
              std::shared_ptr<ROL::OptimizationProblem<double>> opt,
              std::shared_ptr<Teuchos::ParameterList> parlist) {
             new (&instance) ROL::OptimizationSolver<double>(*opt, *parlist);
           })
      .def("solve", [](ROL::OptimizationSolver<double>& instance) {
          instance.solve(std::cout);
      });
}
