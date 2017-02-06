
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>
#include <ROL_BoundConstraint.hpp>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"

#include "EigenVector.h"
#include "CustomLA.h"
#include "PyObjective.h"

PYBIND11_PLUGIN(ROL)
{
  py::module m("ROL", "pybind11 ROL plugin");
  
  py::class_<ROL::Vector<double>, std::shared_ptr<ROL::Vector<double>>>(m, "Vector");

  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>, ROL::Vector<double>, std::shared_ptr<ROL::StdVector<double>>>(m, "StdVector")
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
	.def("__getitem__", [](ROL::StdVector<double> &vec, const py::slice& slice)
	  {
		Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
		py::size_t start, stop, step, slicelength;
		if (!slice.compute((py::size_t)(vvec->size()), &start, &stop, &step, &slicelength))
		  throw py::error_already_set();
		auto res = std::make_shared<std::vector<double>>();
		res->resize(slicelength);
		for (int i = start; i < stop; i=i+step) {
		  res->push_back((*vvec)[i]);
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
  py::class_<EigenVector, ROL::Vector<double>, std::shared_ptr<EigenVector>>(m, "EigenVector", py::buffer_protocol())
	.def(py::init<const int>())
	.def(py::init<const py::array_t<double>>())
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

  py::class_<CustomLA, ROL::Vector<double>, std::shared_ptr<CustomLA>>(m, "CustomLA")
  //py::class_<CustomLA, std::shared_ptr<CustomLA>>(m, "CustomLA")
	.def(py::init<>())
	//.def("dimension", &EigenVector::dimension)
	.def("clone", &CustomLA::clone)
	.def("norm", &CustomLA::norm)
	.def("checkVector", [](std::shared_ptr<CustomLA>& x, std::shared_ptr<CustomLA>& y, std::shared_ptr<CustomLA>& z)
	  {
	    x->checkVector(*y, *z, true, std::cout);
	  })
    .def("cppScale", &CustomLA::scale)
    .def("cppClone", &CustomLA::clone);

  // ROL::Objective<double>
  //
  py::class_<ROL::Objective<double>, PyObjective, std::shared_ptr<ROL::Objective<double>>> objective(m, "Objective");
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
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj)
	  {
	    instance.run(x, obj, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj, ROL::BoundConstraint<double>& bnd)
	  {
	    instance.run(x, obj, bnd, true, std::cout);
	  });

	// ROL::BoundConstraint
	//
	py::class_<ROL::BoundConstraint<double>, std::shared_ptr<ROL::BoundConstraint<double>>>(m, "BoundConstraint")
      .def("__init__",
	    [](ROL::BoundConstraint<double> &instance, std::shared_ptr<ROL::Vector<double>> x_lo,
		   std::shared_ptr<ROL::Vector<double>> x_up, double scale)
		{
		  new (&instance) ROL::BoundConstraint<double>(Teuchos::rcp(x_lo), Teuchos::rcp(x_up), scale);
		});

	// ROL::MoreauYosidaPenalty<double>
	//
	py::class_<ROL::MoreauYosidaPenalty<double>, ROL::Objective<double>, std::shared_ptr<ROL::MoreauYosidaPenalty<double>>>(m, "MoreauYosidaPenalty")
	  .def("__init__",
			  [](ROL::MoreauYosidaPenalty<double> &instance, std::shared_ptr<ROL::Objective<double>> obj,
				 std::shared_ptr<ROL::BoundConstraint<double>> bnd,
				 ROL::Vector<double>& x,
				 const double mu)
			  {
			    new (&instance) ROL::MoreauYosidaPenalty<double>(Teuchos::rcp(obj), Teuchos::rcp(bnd), x, mu);
		      });


  return m.ptr();
}
