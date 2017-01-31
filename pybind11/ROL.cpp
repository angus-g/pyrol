
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>

PYBIND11_PLUGIN(ROL)
{
  py::module m("ROL", "pybind11 ROL plugin");

  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
           Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n));
           new (&instance) ROL::StdVector<double>(tp);
         })
    .def("norm", &ROL::StdVector<double>::norm)
    .def("dimension", &ROL::StdVector<double>::dimension);

  // ROL::Objective<double>
  //
  py::class_<ROL::Objective<double>>(m, "Objective")
    .def("value", &ROL::Objective<double>::value)
    .def("gradient", &ROL::Objective<double>::gradient);

  // ROL::Algorithm<double>
  //
  py::class_<ROL::Algorithm<double>>(m, "Algorithm")
    .def("__init__",
         [](ROL::Algorithm<double> &instance, const std::string& str, const std::map<std::string, std::string>& params)
         {
           Teuchos::ParameterList param_list;
           new (&instance) ROL::Algorithm<double>(str, param_list);
         });

  return m.ptr();
}
