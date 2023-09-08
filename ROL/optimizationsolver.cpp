#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <map>

#include <cereal/types/base_class.hpp>

#include <ROL_OptimizationProblem.hpp>
#include <ROL_OptimizationSolver.hpp>

void init_optimizationsolver(py::module& m) {
  // ROL::OptimizationSolver<double>
  //

  py::class_<ROL::OptimizationSolver<double>,
             std::shared_ptr<ROL::OptimizationSolver<double>>>(
      m, "OptimizationSolver")
      .def(py::init<ROL::OptimizationProblem<double>&, ROL::ParameterList&>())
      .def("solve",
           [](ROL::OptimizationSolver<double>& instance, bool print_output) {
             if(print_output)
               instance.solve(std::cout);
             else
               instance.solve();
           }, py::arg("print_output")=true)
      .def("getAlgorithmState",
           &ROL::OptimizationSolver<double>::getAlgorithmState);
}
