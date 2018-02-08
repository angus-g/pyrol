#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_AugmentedLagrangian.hpp>

void init_augmentedlagrangian(py::module& m) {
    // ROL::AugmentedLagrangian<double>
    //

    py::class_<ROL::AugmentedLagrangian<double>, ROL::Objective<double>,
               std::shared_ptr<ROL::AugmentedLagrangian<double>>>(
        m, "AugmentedLagrangian")
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Constraint<double>>,
                      ROL::Vector<double>&, double, ROL::Vector<double>&,
                      ROL::Vector<double>&, ROL::ParameterList&>());
}
