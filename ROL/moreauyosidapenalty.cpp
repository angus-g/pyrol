#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Bounds.hpp>
#include <ROL_MoreauYosidaPenalty.hpp>

void init_moreauyosidapenalty(py::module& m) {
    // ROL::MoreauYosidaPenalty<double>
    //
    py::class_<ROL::MoreauYosidaPenalty<double>, ROL::Objective<double>,
               std::shared_ptr<ROL::MoreauYosidaPenalty<double>>>(
        m, "MoreauYosidaPenalty")
        .def(py::init<std::shared_ptr<ROL::Objective<double>>,
                      std::shared_ptr<ROL::Bounds<double>>,
                      ROL::Vector<double>&, const double>());
}
