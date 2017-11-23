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
                      std::shared_ptr<ROL::Vector<double>>>())
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::shared_ptr<ROL::Constraint<double>>,
                      std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econ") = (ROL::Constraint<double>*)nullptr,
             py::arg("emul") = (ROL::Vector<double>*)nullptr)
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::vector<std::shared_ptr<ROL::Constraint<double>>>,
                      std::vector<std::shared_ptr<ROL::Vector<double>>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econ") = (std::vector<ROL::Constraint<double>>*)nullptr,
             py::arg("emul") = (std::vector<ROL::Vector<double>>*)nullptr);
}
