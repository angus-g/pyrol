#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "vector.hpp"

#include <ROL_Vector.hpp>
#include "py_shared_ptr.hpp"
PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

void init_vector(py::module& m) {
  py::class_<ROL::Vector<double>, py_shared_ptr<ROL::Vector<double>>,
             PyVector>(m, "Vector")
      .def(py::init<>())
      .def("checkVector", [](ROL::Vector<double>& instance, ROL::Vector<double>& x,
              ROL::Vector<double>& y) {
        return instance.checkVector(x, y);
      });
}
