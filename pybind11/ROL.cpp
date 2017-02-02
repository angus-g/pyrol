
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>

#include "EigenVector.h"

class PyObjective : public ROL::Objective<double>
{
	public:
		//using ROL::Objective<double>::Objective;

		double value( const ROL::Vector<double> &x, double &tol ) override
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "value");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(x, tol).cast<double>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'value'.");
		}
		void gradient( ROL::Vector<double> &g, const ROL::Vector<double> &x, double &tol ) override
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "gradient");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(g, x, tol).cast<void>();
			else
				ROL::Objective<double>::gradient(g, x, tol);
		}
};


PYBIND11_PLUGIN(ROL)
{
  py::module m("ROL", "pybind11 ROL plugin");

  py::class_<ROL::Vector<double>>(m, "Vector");

  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>, ROL::Vector<double>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
           Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n, 0.0));
           new (&instance) ROL::StdVector<double>(tp);
         })
    .def("norm", &ROL::StdVector<double>::norm)
    .def("dimension", &ROL::StdVector<double>::dimension)
	.def("__setitem__", [](ROL::StdVector<double> &vec, const int& idx, const double& val)
	  {
	  Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
		if(idx >= vvec->size())
		{
		  throw py::index_error();
		}else
		{
	      (*vvec)[idx] = val;
		}
	  }
	)
	.def("__getitem__", [](ROL::StdVector<double> &vec, const int& idx)
	  {
	  Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
		if(idx >= vvec->size())
		{
		  throw py::index_error();
		}else
		{
	      return (*vvec)[idx];
		}
	  }
	)
	.def("scale", &ROL::StdVector<double>::scale);

  // EigenVector
  //
  py::class_<EigenVector, ROL::Vector<double>>(m, "EigenVector")
    .def("__init__",
         [](EigenVector &instance, int n) {
           new (&instance) EigenVector(n);
         })
	.def("norm", &EigenVector::norm)
	.def("dimension", &EigenVector::dimension)
	.def("__setitem__", [](EigenVector &vec, const int& idx, const double& val)
	  {
	  auto vvec = vec.getVector();
		if(idx >= vvec->size())
		{
		  throw py::index_error();
		}else
		{
		  (*vvec)[idx] = val;
		}
	  }
	)
	.def("__getitem__", [](EigenVector &vec, const int& idx)
	  {
	  auto vvec = vec.getVector();
		if(idx >= vvec->size())
		{
		  throw py::index_error();
		}else
		{
		  return (*vvec)[idx];
		}
	  }
	)
	.def("scale", &EigenVector::scale);

  // ROL::Objective<double>
  //
  py::class_<ROL::Objective<double>, PyObjective> objective(m, "Objective");
  objective.def(py::init<>())
    .def("value", &ROL::Objective<double>::value)
    .def("gradient", &ROL::Objective<double>::gradient);

  // ROL::Algorithm<double>
  //
  py::class_<ROL::Algorithm<double>>(m, "Algorithm")
    .def("__init__",
         [](ROL::Algorithm<double> &instance, const std::string& str, const std::string& xml_params)
         {
           Teuchos::RCP<Teuchos::ParameterList> params =
             Teuchos::getParametersFromXmlString(xml_params);

           new (&instance) ROL::Algorithm<double>(str, *params);
         })
    //.def("run", (std::vector<std::string> (ROL::Algorithm<double>::*)(ROL::Vector<double>&, ROL::Objective<double>&, bool, std::ostream&)) &ROL::Algorithm<double>::run);
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj)
	  {
	    instance.run(x, obj, true, std::cout);
	  });

  return m.ptr();
}
