#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>

#include <ROL_NewtonKrylov_U.hpp>

void init_descentdirection(py::module& m) {
  py::class_<ROL::DescentDirection_U<double>, std::shared_ptr<ROL::DescentDirection_U<double>>>(m, "DescentDirection");

  py::class_<ROL::NewtonKrylov_U<double>, ROL::DescentDirection_U<double>, std::shared_ptr<ROL::NewtonKrylov_U<double>>>(m, "NewtonKrylov")
    .def(py::init<ROL::ParameterList&>());
}
