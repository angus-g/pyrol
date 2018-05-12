#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Bounds.hpp>
#include <iostream>

#include "py_shared_ptr.hpp"
PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

class PyBounds : public ROL::Bounds<double> {
   public:
    // PyBounds()
    //    : ROL::Bounds<double>()
    //{}

    PyBounds(const std::shared_ptr<ROL::Vector<double>>& x_lo,
             const std::shared_ptr<ROL::Vector<double>>& x_up,
             const double scale = 1)
        : ROL::Bounds<double>(x_lo, x_up, scale) {}

    PyBounds(const ROL::Vector<double>& x,
             bool isLower = true,
             double scale = 1)
        : ROL::Bounds<double>(x, isLower, scale) {}

    void project(ROL::Vector<double>& x) override {
        PYBIND11_OVERLOAD(void, ROL::Bounds<double>, project, x);
    }
};

void init_bounds(py::module& m) {
    py::class_<ROL::BoundConstraint<double>,
               py_shared_ptr<ROL::BoundConstraint<double>>>(
        m, "BoundConstraint");

    //
    // ROL::Bounds
    //
    py::class_<ROL::Bounds<double>, ROL::BoundConstraint<double>, PyBounds,
               std::shared_ptr<ROL::Bounds<double>>>(m, "Bounds")
        .def(py::init<const std::shared_ptr<ROL::Vector<double>>&,
                      const std::shared_ptr<ROL::Vector<double>>&, double>(), 
                      py::arg("x_lo"), py::arg("x_up"), py::arg("scale") = 1.0)
        .def(py::init<const ROL::Vector<double>&, bool, double>(),
            py::arg("x_lo"), py::arg("isLower")=true, py::arg("scale") = 1.0)
        .def("test", [](const ROL::Bounds<double> &inst){
            std::cout << "lower dim:" << inst.getLowerBound()->dimension() << std::endl;
            std::cout << "upper dim:" << inst.getUpperBound()->dimension() << std::endl;
            });
}
