#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_MoreauYosidaPenalty.hpp>
#include <ROL_Bounds.hpp>

void init_moreauyosidapenalty(py::module& m) {
    // ROL::MoreauYosidaPenalty<double>
    //
    py::class_<ROL::MoreauYosidaPenalty<double>, ROL::Objective<double>,
               std::shared_ptr<ROL::MoreauYosidaPenalty<double>>>(
        m, "MoreauYosidaPenalty")
        .def("__init__", [](ROL::MoreauYosidaPenalty<double>& instance,
                            std::shared_ptr<ROL::Objective<double>> obj,
                            std::shared_ptr<ROL::Bounds<double>> bnd,
                            ROL::Vector<double>& x, const double mu) {
            new (&instance) ROL::MoreauYosidaPenalty<double>(obj, bnd, x, mu);
        });
}
