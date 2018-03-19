#include <pybind11/iostream.h>
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
             py::scoped_ostream_redirect stream(
                 std::cout,                                // std::ostream&
                 py::module::import("sys").attr("stdout")  // Python output
                 );
             instance.solve(std::cout);
           })
      .def("getAlgorithmState",
           &ROL::OptimizationSolver<double>::getAlgorithmState);
}
