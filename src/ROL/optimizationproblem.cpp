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
                      std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("obj"), py::arg("sol"))
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>>(),
             py::arg("obj"), py::arg("sol"), py::arg("bnd"))
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::shared_ptr<ROL::Constraint<double>>,
                      std::shared_ptr<ROL::Vector<double>>>(),
             py::arg("obj"), py::arg("sol"), py::arg("bnd"), py::arg("econ"),
             py::arg("emul"))
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Vector<double>>,
                      std::shared_ptr<ROL::BoundConstraint<double>>,
                      std::vector<std::shared_ptr<ROL::Constraint<double>>>,
                      std::vector<std::shared_ptr<ROL::Vector<double>>>>(),
             py::arg("obj"), py::arg("sol"), py::arg("bnd"), py::arg("econ"),
             py::arg("emul"));
}
