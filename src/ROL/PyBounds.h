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
        PYBIND11_OVERLOAD(void, ROL::Bounds<double>, project, x);
    }
};
