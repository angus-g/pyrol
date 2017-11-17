
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <stdexcept>
namespace py = pybind11;

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Types.hpp>
#include "Teuchos_ParameterList.hpp"

#include "CustomLA.h"
#include "PyObjective.h"
#include "PyConstraint.h"
#include "PyBounds.h"

void dictToParameterList(py::dict param_dict, Teuchos::ParameterList &parlist)
{
    for(auto item : param_dict)
    {
        auto key = item.first;
        auto value = item.second;
        if(py::str(key, true).check())
        {
            if(py::str(value, true).check())
                parlist.set(std::string(py::str(key)), std::string(py::str(value)));
            else if(py::bool_(value, true).check())
                parlist.set(std::string(py::str(key)), value.cast<bool>());
            else if(py::int_(value, true).check())
                parlist.set(std::string(py::str(key)), value.cast<int>());
            else if(py::float_(value, true).check())
                parlist.set(std::string(py::str(key)), value.cast<double>());
            else if(py::dict(value, true).check())
            {
                auto &sublist = parlist.sublist(std::string(py::str(key)));
                dictToParameterList(value.cast<py::dict>(), sublist);
            }
            else
                throw std::runtime_error(std::string("Don't know what to do with value."));
        }
        else
            throw std::runtime_error(std::string("Don't know what to do with key."));
    }
}
PYBIND11_PLUGIN(_ROL)
{
  py::module m("_ROL", "PyROL provides Python wrappers for a subset of the Trilinos ROL library.");
  m.attr("__version__") = "0.1.1";

  py::class_<ROL::Vector<double>, std::shared_ptr<ROL::Vector<double>>>(m, "Vector");

  //
  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>, ROL::Vector<double>, std::shared_ptr<ROL::StdVector<double>>>(m, "StdVector")
    .def("__init__",
         [](ROL::StdVector<double> &instance, int n) {
            auto tp = std::make_shared<std::vector<double>>(n, 0.0);
           new (&instance) ROL::StdVector<double>(tp);
         })
    .def("norm", &ROL::StdVector<double>::norm, "L2 norm of the vector")
    .def("dimension", &ROL::StdVector<double>::dimension, "Size of the vector")
	.def("__setitem__", [](ROL::StdVector<double> &vec, const int& idx, const double& val)
	  {
	  auto vvec = vec.getVector();
          if(idx >= (int)vvec->size())
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
		auto vvec = vec.getVector();
		py::size_t start, stop, step, slicelength;
		if (!slice.compute((py::size_t)(vvec->size()), &start, &stop, &step, &slicelength))
		  throw py::error_already_set();
		auto res = std::make_shared<std::vector<double>>();
		res->resize(slicelength);
		for (unsigned int i = start; i < stop; i=i+step) {
		  res->push_back((*vvec)[i]);
		}
	  }
	)
	.def("__getitem__", [](ROL::StdVector<double> &vec, const int& idx)
	  {
	  auto vvec = vec.getVector();
          if(idx >= (int)vvec->size())
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
	//.def("clone", &CustomLA::clone)
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
  py::class_<ROL::Objective<double>, PyObjective, std::shared_ptr<ROL::Objective<double>>> objective(m, "Objective", "Base class for the objective class. Python objectives need to inherit from this class.");
  objective.def(py::init<>())
    //.def("value", &ROL::Objective<double>::value)
    //.def("gradient", &ROL::Objective<double>::gradient)
    //.def("update", &ROL::Objective<double>::update)
    .def("checkGradient", [](ROL::Objective<double>& instance, ROL::Vector<double>& x, ROL::Vector<double>& d, int steps, int order)->std::vector<double>
	  {
	    auto res = instance.checkGradient(x, d, true, std::cout, steps, order);
        return res[0];
	  });

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
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj, ROL::Bounds<double>& bnd)
	  {
	    instance.run(x, obj, bnd, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Vector<double>& l, ROL::Objective<double>& obj, ROL::Constraint<double>& con)
	  {
	    instance.run(x, l, obj, con, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Vector<double>& l, ROL::Objective<double>& obj, ROL::Constraint<double>& con, ROL::Bounds<double>& bnd)
	  {
	    instance.run(x, l, obj, con, bnd, true, std::cout);
	  })
    .def("run", [](ROL::Algorithm<double> &instance, ROL::OptimizationProblem<double>& opt) {
	    instance.run(opt, true, std::cout);
	  })
    
    .def("get_state", [](ROL::Algorithm<double> &instance)
    {
        return instance.getState();   
    });

    py::class_<ROL::AlgorithmState<double>, std::shared_ptr<ROL::AlgorithmState<double>>>(m, "AlgorithmState")
        .def_readwrite("gnorm", &ROL::AlgorithmState<double>::gnorm)
        .def_readwrite("cnorm", &ROL::AlgorithmState<double>::cnorm)
        .def_readwrite("snorm", &ROL::AlgorithmState<double>::snorm)
        .def_readwrite("value", &ROL::AlgorithmState<double>::value);

	// ROL::Bounds
	//
    py::class_<ROL::Bounds<double>, PyBounds, std::shared_ptr<ROL::Bounds<double>>>(m, "Bounds")
	//py::class_<ROL::Bounds<double>, std::shared_ptr<ROL::Bounds<double>>>(m, "Bounds")
      //.def(py::init<>())
      .def("__init__",
	    [](ROL::Bounds<double> &instance, std::shared_ptr<ROL::Vector<double>> x_lo,
		   std::shared_ptr<ROL::Vector<double>> x_up, double scale)
		{
		  //new (&instance) ROL::Bounds<double>(Teuchos::rcp(x_lo), Teuchos::rcp(x_up), scale);
		  new (&instance) PyBounds(x_lo, x_up, scale);
		})
      .def("test", [](ROL::Bounds<double> &instance, ROL::Vector<double>& x)
        {
          instance.project(x);
        });
      //.def("__init__",
      //  [](ROL::Bounds<double> &instance)
      //  {
      //    new (&instance) PyBounds();
      //    instance.deactivate();
      //  });

	// ROL::MoreauYosidaPenalty<double>
	//
	py::class_<ROL::MoreauYosidaPenalty<double>, ROL::Objective<double>, std::shared_ptr<ROL::MoreauYosidaPenalty<double>>>(m, "MoreauYosidaPenalty")
	  .def("__init__",
			  [](ROL::MoreauYosidaPenalty<double> &instance, std::shared_ptr<ROL::Objective<double>> obj,
				 std::shared_ptr<ROL::Bounds<double>> bnd,
				 ROL::Vector<double>& x,
				 const double mu)
			  {
			    new (&instance) ROL::MoreauYosidaPenalty<double>(obj, bnd, x, mu);
		      });

    // ROL::Constraint<double>
	//
	py::class_<ROL::Constraint<double>, PyConstraint, std::shared_ptr<ROL::Constraint<double>>>(m, "Constraint")
	  .def(py::init<>())
	  //.def("value", &ROL::Constraint<double>::value)
	  //.def("applyJacobian", &ROL::Constraint<double>::applyJacobian)
	  //.def("applyAdjointJacobian", [](ROL::Constraint<double>& instance,
	  //  ROL::Vector<double> &ajv, const ROL::Vector<double> &v,
	  //  const ROL::Vector<double> &x, double &tol)
	  //  {
	  //    instance.applyAdjointJacobian(ajv, v, x, tol);
	  //  })
	  .def("checkApplyJacobian", [](ROL::Constraint<double>& instance,
        ROL::Vector<double>& x, ROL::Vector<double>& v, ROL::Vector<double>& jv,
        int steps, int order)
        {
          instance.checkApplyJacobian(x, v, jv, true, std::cout, steps, order);
        })
      .def("checkAdjointConsistencyJacobian", [](ROL::Constraint<double>& instance,
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
		  std::shared_ptr<ROL::Constraint<double>> con,
		  ROL::Vector<double>& multiplier,
		  double penaltyParameter,
		  ROL::Vector<double>& optVec,
		  ROL::Vector<double>& conVec,
		  Teuchos::ParameterList& params)
		{
		  new (&instance) ROL::AugmentedLagrangian<double>(obj,
				  con,
				  multiplier,
				  penaltyParameter,
				  optVec,
				  conVec,
				  params);
		});

    // ROL::OptimizationProblem<double>
    //

    //py::class_<ROL::OptimizationProblem<double>, std::shared_ptr<ROL::OptimizationProblem<double>>>(m, "OptimizationProblem")
    //  .def("__init__",
    //    [](ROL::OptimizationProblem<double>& instance,
    //       std::shared_ptr<ROL::Objective<double>> obj,
    //       std::shared_ptr<ROL::Vector<double>> sol,
    //       std::shared_ptr<ROL::Bounds<double>> bnd,
    //       std::shared_ptr<Teuchos::ParameterList> params)
    //    {
    //      new (&instance) ROL::OptimizationProblem<double>(obj,
    //              sol,
    //              bnd,
	//              params);
	//    });

    py::class_<Teuchos::ParameterList, std::shared_ptr<Teuchos::ParameterList>>(m, "ParameterList",
                                                   "Create a ParameterList object from an XML string")
        .def("__init__",
             [](Teuchos::ParameterList& instance,
                std::string xml_params)
             {
               new (&instance) Teuchos::ParameterList(*(Teuchos::getParametersFromXmlString(xml_params)));
             })
        .def("__init__",
                [](Teuchos::ParameterList& instance,
                   py::dict param_dict,
                   std::string name)
            {
                auto res = new (&instance) Teuchos::ParameterList(name);
                dictToParameterList(param_dict, *res);
            })
        .def("print", [](Teuchos::ParameterList &instance){instance.print();});


  return m.ptr();
}
