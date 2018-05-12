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
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::Constraint<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econ") = (ROL::Constraint<double>*)nullptr,
             py::arg("emul") = (ROL::Vector<double>*)nullptr,
             py::arg("icon") = (ROL::Constraint<double>*)nullptr,
             py::arg("imul") = (ROL::Vector<double>*)nullptr,
             py::arg("ibnd") = (ROL::BoundConstraint<double>*)nullptr)
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::vector<std::shared_ptr<ROL::Constraint<double>>>,
                      std::vector<std::shared_ptr<ROL::Vector<double>>>,
                      std::vector<std::shared_ptr<ROL::Constraint<double>>>,
                      std::vector<std::shared_ptr<ROL::Vector<double>>>,
                      std::vector<std::shared_ptr<ROL::BoundConstraint<double>>>>(),
             py::arg("obj"), py::arg("sol"),
             py::arg("bnd") = (ROL::BoundConstraint<double>*)nullptr,
             py::arg("econs") = std::vector<std::shared_ptr<ROL::Constraint<double>>>(),
             py::arg("emuls") = std::vector<std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("icons") = std::vector<std::shared_ptr<ROL::Constraint<double>>>(),
             py::arg("imuls") = std::vector<std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("ibnds") = std::vector<std::shared_ptr<ROL::BoundConstraint<double>>>());
}
