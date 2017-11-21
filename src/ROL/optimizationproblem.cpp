#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_OptimizationProblem.hpp>

void init_optimizationproblem(py::module& m) {
    // ROL::OptimizationProblem<double>
    //

    py::class_<ROL::OptimizationProblem<double>,
               std::shared_ptr<ROL::OptimizationProblem<double>>>(
        m, "OptimizationProblem")
        .def("__init__",
             [](ROL::OptimizationProblem<double>& instance,
                std::shared_ptr<ROL::Objective<double>> obj,
                std::shared_ptr<ROL::Vector<double>> sol, py::kwargs kwargs) {

                 std::shared_ptr<ROL::BoundConstraint<double>> bnd = nullptr;
                 std::shared_ptr<ROL::Constraint<double>> econ = nullptr;
                 std::shared_ptr<ROL::Vector<double>> emul = nullptr;

                 if (kwargs) {
                     for (auto item : kwargs) {
                         auto key = item.first.cast<std::string>();
                         if (key == "bnd") {
                             bnd = item.second.cast<std::shared_ptr<
                                 ROL::BoundConstraint<double>>>();
                         } else if (key == "econ") {
                             econ = item.second.cast<
                                 std::shared_ptr<ROL::Constraint<double>>>();
                         } else if (key == "emul") {
                             emul = item.second.cast<
                                 std::shared_ptr<ROL::Vector<double>>>();
                         }
                     }
                 }
                 new (&instance) ROL::OptimizationProblem<double>(obj, sol, bnd,
                                                                  econ, emul);
             },
             "Creates an OptimizationProblem object. \n"
             "Arguments:\nobj: Objective\nsol: Vector\nOptional "
             "Arguments:\nbnd: "
             "BoundConstraint\necon: Constraint\nemul: Vector");
}
