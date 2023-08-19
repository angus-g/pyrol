#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

#include "vector.hpp"
#include <ROL_Types.hpp>

#include <ROL_Secant.hpp>

CEREAL_REGISTER_TYPE(ROL::lBFGS<double>);

void init_secant(py::module& m) {
  py::class_<ROL::Secant<double>, std::shared_ptr<ROL::Secant<double>>>(m, "Secant");

  py::class_<ROL::lBFGS<double>, ROL::Secant<double>, std::shared_ptr<ROL::lBFGS<double>>>(m, "lBFGS")
    .def(py::init<>())
    .def(py::init<int>());
}
