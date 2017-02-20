
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>
#include <ROL_BoundConstraint.hpp>
#include <ROL_EqualityConstraint.hpp>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CustomLA.h"
#include "PyObjective.h"
#include "PyEqualityConstraint.h"

PYBIND11_PLUGIN(ROL)
{
  py::module m("ROL", "PyROL provides Python wrappers for a subset of the Trilinos ROL library.");

  py::class_<ROL::Vector<double>, std::shared_ptr<ROL::Vector<double>>>(m, "Vector");

  //
  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>, ROL::Vector<double>, std::shared_ptr<ROL::StdVector<double>>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
           Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n, 0.0));
           new (&instance) ROL::StdVector<double>(tp);
         })
    .def("norm", &ROL::StdVector<double>::norm, "L2 norm of the vector")
    .def("dimension", &ROL::StdVector<double>::dimension, "Size of the vector")
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
    .def("scale", &ROL::StdVector<double>::scale, "Multiply the vector by a scalar")
        .def("checkVector", [](std::shared_ptr<ROL::StdVector<double>>& x,
                               std::shared_ptr<ROL::StdVector<double>>& y,
                               std::shared_ptr<ROL::StdVector<double>>& z)->std::vector<double>
             {
               return x->checkVector(*y, *z, true, std::cout);
             }, "Check the accuracy of the linear algebra implementation");

  //
  // ROL::CustomLA
  //
  py::class_<CustomLA, ROL::Vector<double>, std::shared_ptr<CustomLA>>(m, "CustomLA", "Custom Vector and Linear Algebra base class")
  //py::class_<CustomLA, std::shared_ptr<CustomLA>>(m, "CustomLA")
	.def(py::init<>())
	.def("clone", &CustomLA::clone)
	.def("norm", &CustomLA::norm)
        .def("checkVector", [](std::shared_ptr<CustomLA>& x, std::shared_ptr<CustomLA>& y, std::shared_ptr<CustomLA>& z)->std::vector<double>
	  {
	    return x->checkVector(*y, *z, true, std::cout);
	  });
  //    .def("cppScale", &CustomLA::scale)
  //    .def("cppClone", &CustomLA::clone);

  //
  // ROL::Objective<double>
  //
  py::class_<ROL::Objective<double>, PyObjective, std::shared_ptr<ROL::Objective<double>>> objective(m, "Objective");
  objective.def(py::init<>())
    .def("value", &ROL::Objective<double>::value)
    .def("gradient", &ROL::Objective<double>::gradient)
    .def("update", &ROL::Objective<double>::update);

  // ROL::Algorithm<double>
  //
  py::class_<ROL::Algorithm<double>>(m, "Algorithm")
    .def("__init__",
      [](ROL::Algorithm<double> &instance, const std::string& str, Teuchos::ParameterList& params)
      {
        new (&instance) ROL::Algorithm<double>(str, params);
      })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj)
	  {
	    instance.run(x, obj, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj, ROL::BoundConstraint<double>& bnd)
	  {
	    instance.run(x, obj, bnd, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Vector<double>& l, ROL::Objective<double>& obj, ROL::EqualityConstraint<double>& con)
	  {
	    instance.run(x, l, obj, con, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Vector<double>& l, ROL::Objective<double>& obj, ROL::EqualityConstraint<double>& con, ROL::BoundConstraint<double>& bnd)
	  {
	    instance.run(x, l, obj, con, bnd, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::OptimizationProblem<double>& opt) {
	    instance.run(opt, true, std::cout);
	  });

	// ROL::BoundConstraint
	//
	py::class_<ROL::BoundConstraint<double>, std::shared_ptr<ROL::BoundConstraint<double>>>(m, "BoundConstraint")
      .def("__init__",
	    [](ROL::BoundConstraint<double> &instance, std::shared_ptr<ROL::Vector<double>> x_lo,
		   std::shared_ptr<ROL::Vector<double>> x_up, double scale)
		{
		  new (&instance) ROL::BoundConstraint<double>(Teuchos::rcp(x_lo), Teuchos::rcp(x_up), scale);
		})
      .def("__init__",
	    [](ROL::BoundConstraint<double> &instance)
		{
		  new (&instance) ROL::BoundConstraint<double>();
		  instance.deactivate();
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

    // ROL::EqualityConstraint<double>
	//
	py::class_<ROL::EqualityConstraint<double>, PyEqualityConstraint, std::shared_ptr<ROL::EqualityConstraint<double>>>(m, "EqualityConstraint")
	  .def(py::init<>())
	  .def("value", &ROL::EqualityConstraint<double>::value)
	  .def("applyJacobian", &ROL::EqualityConstraint<double>::applyJacobian)
	  .def("applyAdjointJacobian", [](ROL::EqualityConstraint<double>& instance,
	    ROL::Vector<double> &ajv, const ROL::Vector<double> &v,
		const ROL::Vector<double> &x, double &tol)
		{
		  instance.applyAdjointJacobian(ajv, v, x, tol);
		})
	  .def("checkApplyJacobian", [](ROL::EqualityConstraint<double>& instance,
        ROL::Vector<double>& x, ROL::Vector<double>& v, ROL::Vector<double>& jv,
        int steps, int order)
        {
          instance.checkApplyJacobian(x, v, jv, true, std::cout, steps, order);
        })
      .def("checkAdjointConsistencyJacobian", [](ROL::EqualityConstraint<double>& instance,
        ROL::Vector<double>& w, ROL::Vector<double>& v, ROL::Vector<double>& x)
        {
          instance.checkAdjointConsistencyJacobian(w, v, x, true, std::cout);
        });

	// ROL::AugmentedLagrangian<double>
	//

	py::class_<ROL::AugmentedLagrangian<double>, ROL::Objective<double>, std::shared_ptr<ROL::AugmentedLagrangian<double>>>(m, "AugmentedLagrangian")
	  .def("__init__",
	    [](ROL::AugmentedLagrangian<double>& instance,
		  std::shared_ptr<ROL::Objective<double>> obj,
		  std::shared_ptr<ROL::EqualityConstraint<double>> con,
		  ROL::Vector<double>& multiplier,
		  double penaltyParameter,
		  ROL::Vector<double>& optVec,
		  ROL::Vector<double>& conVec,
		  Teuchos::ParameterList& params)
		{
		  new (&instance) ROL::AugmentedLagrangian<double>(Teuchos::rcp(obj),
				  Teuchos::rcp(con),
				  multiplier,
				  penaltyParameter,
				  optVec,
				  conVec,
				  params);
		});

    // ROL::OptimizationProblem<double>
    //

    py::class_<ROL::OptimizationProblem<double>, std::shared_ptr<ROL::OptimizationProblem<double>>>(m, "OptimizationProblem")
      .def("__init__",
        [](ROL::OptimizationProblem<double>& instance,
           std::shared_ptr<ROL::Objective<double>> obj,
           std::shared_ptr<ROL::Vector<double>> sol,
           std::shared_ptr<ROL::BoundConstraint<double>> bnd,
           std::shared_ptr<Teuchos::ParameterList> params)
        {
          new (&instance) ROL::OptimizationProblem<double>(Teuchos::rcp(obj),
                  Teuchos::rcp(sol),
                  Teuchos::rcp(bnd),
				  Teuchos::rcp(params));
		});

    py::class_<Teuchos::ParameterList, std::shared_ptr<Teuchos::ParameterList>>(m, "ParameterList",
                                                   "Create a ParameterList object from an XML string")
        .def("__init__",
             [](Teuchos::ParameterList& instance,
                std::string xml_params)
             {
               new (&instance) Teuchos::ParameterList(*(Teuchos::getParametersFromXmlString(xml_params)));
             });

  return m.ptr();
}
