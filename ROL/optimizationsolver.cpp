#include <pybind11/pybind11.h>

#include <ROL_OptimizationProblem.hpp>
#include <ROL_OptimizationSolver.hpp>

namespace py = pybind11;
void init_optimizationsolver(py::module& m) {
  // ROL::OptimizationSolver<double>
  //

  py::class_<ROL::OptimizationSolver<double>,
             std::shared_ptr<ROL::OptimizationSolver<double>>>(
      m, "OptimizationSolver")
      .def(py::init<ROL::OptimizationProblem<double>&, ROL::ParameterList&>())
      .def("solve",
           [](ROL::OptimizationSolver<double>& instance) {
             instance.solve(std::cout);
           }, 
      py::call_guard<py::scoped_ostream_redirect,
                     py::scoped_estream_redirect>())
      .def("getAlgorithmState",
           &ROL::OptimizationSolver<double>::getAlgorithmState);
}
