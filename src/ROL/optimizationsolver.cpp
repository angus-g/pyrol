#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

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
        .def(py::init<ROL::OptimizationProblem<double>&,
                      Teuchos::ParameterList&>())
        .def("solve", [](ROL::OptimizationSolver<double>& instance) {
            instance.solve(std::cout);
        });
}
