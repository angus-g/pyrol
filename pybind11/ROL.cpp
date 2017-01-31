
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>

PYBIND11_PLUGIN(ROL)
{
  py::module m("ROL", "pybind11 ROL plugin");

  py::class_<ROL::Vector<double>>(m, "Vector");

  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
           Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n, 1.0));
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
           Teuchos::ParameterList parlist;
           parlist.sublist("Step").sublist("Line Search")
             .sublist("Descent Method").set("Type", "Newton-Krylov");
           parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
           parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
           parlist.sublist("Status Test").set("Iteration Limit",100);

           new (&instance) ROL::Algorithm<double>(str, parlist);
         })
    .def("run", (std::vector<std::string> (ROL::Algorithm<double>::*)(ROL::Vector<double>&, ROL::Objective<double>&, bool, std::ostream&)) &ROL::Algorithm<double>::run);

  return m.ptr();
}
