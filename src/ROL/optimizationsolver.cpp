#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_OptimizationSolver.hpp>

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
                 new (&instance)
                     ROL::OptimizationSolver<double>(*opt, *parlist);
             })
        .def("solve", [](ROL::OptimizationSolver<double>& instance) {
            instance.solve(std::cout);
        });
}
