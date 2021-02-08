#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_ParameterList.hpp>
#include <ROL_StatusTest.hpp>

void init_statustest(py::module& m) {
  py::class_<ROL::StatusTest<double>, std::shared_ptr<ROL::StatusTest<double>>>(m, "StatusTest")
    .def(py::init<ROL::ParameterList&>());
}
