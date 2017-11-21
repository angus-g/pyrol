#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_AugmentedLagrangian.hpp>

void init_augmentedlagrangian(py::module& m) {
    // ROL::AugmentedLagrangian<double>
    //

    py::class_<ROL::AugmentedLagrangian<double>, ROL::Objective<double>,
               std::shared_ptr<ROL::AugmentedLagrangian<double>>>(
        m, "AugmentedLagrangian")
        .def("__init__",
             [](ROL::AugmentedLagrangian<double>& instance,
                std::shared_ptr<ROL::Objective<double>> obj,
                std::shared_ptr<ROL::Constraint<double>> con,
                ROL::Vector<double>& multiplier, double penaltyParameter,
                ROL::Vector<double>& optVec, ROL::Vector<double>& conVec,
                Teuchos::ParameterList& params) {
                 new (&instance) ROL::AugmentedLagrangian<double>(
                     obj, con, multiplier, penaltyParameter, optVec, conVec,
                     params);
             });
}
