#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <ROL_OptimizationProblem.hpp>

void init_optimizationproblem(py::module& m) {
    // ROL::OptimizationProblem<double>
    //

    py::class_<ROL::OptimizationProblem<double>,
               std::shared_ptr<ROL::OptimizationProblem<double>>>(
        m, "OptimizationProblem")
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::shared_ptr<ROL::Constraint<double>>,
                      std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econ") = (ROL::Constraint<double>*)nullptr,
             py::arg("emul") = (ROL::Vector<double>*)nullptr)
        // Have to call the parameters econs and emuls when accepting lists.
        // Otherwise the wrong function is called and we get an error.
        // Should file a bug report for this but haven't been able to reduce it
        // to a MFE.
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::vector<std::shared_ptr<ROL::Constraint<double>>>,
                      std::vector<std::shared_ptr<ROL::Vector<double>>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econs") = std::vector<std::shared_ptr<ROL::Constraint<double>>>(),
             py::arg("emuls") = std::vector<std::shared_ptr<ROL::Vector<double>>>());
}
