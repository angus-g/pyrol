
#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_StdVector.hpp>

PYBIND11_PLUGIN(ROL)
{

  py::module m("ROL", "pybind11 ROL plugin");

  py::class_<ROL::StdVector<double>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
           Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n));
           new (&instance) ROL::StdVector<double>(tp);
         })
    .def("dimension", &ROL::StdVector<double>::dimension);


  return m.ptr();

}
