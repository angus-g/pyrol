#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "py_shared_ptr.hpp"
PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

#include <ROL_StatusTest.hpp>

class PyStatusTest : public ROL::StatusTest<double> {
public:
  using ROL::StatusTest<double>::StatusTest;

  bool check(ROL::AlgorithmState<double> &state) {
    PYBIND11_OVERRIDE(bool, ROL::StatusTest<double>, check, state);
  }
};

void init_statustest(py::module &m) {
  py::class_<ROL::StatusTest<double>, PyStatusTest, py_shared_ptr<ROL::StatusTest<double>>>(m, "StatusTest")
    .def(py::init<ROL::ParameterList&>())
    .def("check", &ROL::StatusTest<double>::check);
}
