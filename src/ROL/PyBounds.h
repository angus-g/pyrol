#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Bounds.hpp>

class PyBounds : public ROL::Bounds<double> {
 public:
  // PyBounds()
  //    : ROL::Bounds<double>()
  //{}

  PyBounds(const std::shared_ptr<ROL::Vector<double> > &x_lo,
           const std::shared_ptr<ROL::Vector<double> > &x_up,
           const double scale = 1)
      : ROL::Bounds<double>(x_lo, x_up, scale) {}

  void project(ROL::Vector<double> &x) override {
    py::gil_scoped_acquire gil;
    py::function overload = py::get_overload(this, "project");
    if (overload)
      return overload.operator()<py::return_value_policy::reference>(x)
          .cast<void>();
    else
      ROL::Bounds<double>::project(x);
  }
};
