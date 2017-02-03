
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <Python.h>
namespace py = pybind11;

#include <ROL_StdVector.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Algorithm.hpp>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"

#include "EigenVector.h"
class CustomLA : public std::enable_shared_from_this<CustomLA>, public ROL::Vector<double>
{
	public:
		CustomLA()
		{
		}

		virtual void plus(const ROL::Vector<double>& x)
		{
			const CustomLA& xx = Teuchos::dyn_cast<const CustomLA>(x);
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "plus");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(xx).cast<void>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'plus'.");
		}

		virtual void scale(const double alpha)
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "scale");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(alpha).cast<void>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'scale'.");
		}

		virtual double dot(const ROL::Vector<double> &x) const
		{
			const CustomLA& xx = Teuchos::dyn_cast<const CustomLA>(x);
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "dot");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(xx).cast<double>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'dot'.");
		}

		virtual double norm() const
		{
			return std::sqrt(this->dot(*this));
		}

		//virtual int dimension() const {
		//    return _vec->size();
		//}

		virtual Teuchos::RCP<ROL::Vector<double>> clone() const
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "clone");
			if (overload)
			{
				auto res = overload.operator()<py::return_value_policy::reference>();
				Py_INCREF(res.ptr());
				return Teuchos::rcp(res.cast<std::shared_ptr<CustomLA>>());
			}
			else
				py::pybind11_fail("Tried to call pure virtual function 'clone'.");
		}

};

class PyObjective : public ROL::Objective<double>
{
	public:

		double value( const ROL::Vector<double>& x, double &tol ) override
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

  py::class_<ROL::Vector<double>, std::shared_ptr<ROL::Vector<double>>>(m, "Vector");

  // ROL::StdVector<double>
  //
  //py::class_<ROL::StdVector<double>, ROL::Vector<double>>(m, "StdVector")
  //  .def("__init__",
  //       [](ROL::StdVector<double> &instance, int n) {
  //         Teuchos::RCP<std::vector<double>> tp = Teuchos::rcp<std::vector<double>>(new std::vector<double>(n, 0.0));
  //         new (&instance) ROL::StdVector<double>(tp);
  //       })
  //  .def("norm", &ROL::StdVector<double>::norm)
  //  .def("dimension", &ROL::StdVector<double>::dimension)
  //  .def("__setitem__", [](ROL::StdVector<double> &vec, const int& idx, const double& val)
  //    {
  //    Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
  //      if(idx >= vvec->size())
  //      {
  //        throw py::index_error();
  //      }else
  //      {
  //        (*vvec)[idx] = val;
  //      }
  //    }
  //  )
  //  .def("__getitem__", [](ROL::StdVector<double> &vec, const py::slice& slice)
  //    {
  //      Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
  //      py::size_t start, stop, step, slicelength;
  //      if (!slice.compute((py::size_t)(vvec->size()), &start, &stop, &step, &slicelength))
  //        throw py::error_already_set();
  //      auto res = std::make_shared<std::vector<double>>();
  //      res->resize(slicelength);
  //      for (int i = start; i < stop; i=i+step) {
  //        res->push_back((*vvec)[i]);
  //      }
  //    }
  //  )
  //  .def("__getitem__", [](ROL::StdVector<double> &vec, const int& idx)
  //    {
  //    Teuchos::RCP<std::vector<double>> vvec = vec.getVector();
  //      if(idx >= vvec->size())
  //      {
  //        throw py::index_error();
  //      }else
  //      {
  //        return (*vvec)[idx];
  //      }
  //    }
  //  )
  //  .def("scale", &ROL::StdVector<double>::scale);

  // EigenVector
  //
  //py::class_<EigenVector, ROL::Vector<double>>(m, "EigenVector", py::buffer_protocol())
  //  .def(py::init<const int>())
  //  .def(py::init<const py::array_t<double>>())
  //  .def("norm", &EigenVector::norm)
  //  .def("dimension", &EigenVector::dimension)
  //  .def("__setitem__", [](EigenVector &vec, const int& idx, const double& val)
  //    {
  //    auto vvec = vec.getVector();
  //      if(idx >= vvec->size())
  //      {
  //        throw py::index_error();
  //      }else
  //      {
  //        (*vvec)[idx] = val;
  //      }
  //    }
  //  )
  //  .def("__getitem__", [](EigenVector &vec, const int& idx)
  //    {
  //    auto vvec = vec.getVector();
  //      if(idx >= vvec->size())
  //      {
  //        throw py::index_error();
  //      }else
  //      {
  //        return (*vvec)[idx];
  //      }
  //    }
  //  )
  //  .def("scale", &EigenVector::scale);

  //py::class_<CustomLA, ROL::Vector<double>>(m, "CustomLA")
  py::class_<CustomLA, ROL::Vector<double>, std::shared_ptr<CustomLA>>(m, "CustomLA")
	//.def(py::init<std::function<void(CustomLA, CustomLA)>, std::function<void(CustomLA, double)>>())
	.def(py::init<>())
	//.def("norm", &EigenVector::norm)
	//.def("dimension", &EigenVector::dimension)
	//.def("plus", &CustomLA::plus)
	//.def("scale", &CustomLA::scale)
	.def("clone", &CustomLA::clone)
	.def("norm", &CustomLA::norm)
	.def("checkVector", [](std::shared_ptr<CustomLA>& x, std::shared_ptr<CustomLA>& y, std::shared_ptr<CustomLA>& z)
	  {
	    x->checkVector(*y, *z, true, std::cout);
	  })
    .def("checkScale", [](std::shared_ptr<CustomLA>& instance) 
	  {
	    instance->scale(0.5);
	  })
    .def("checkClone", [](std::shared_ptr<CustomLA>& instance) 
	  {
	    instance->scale(0.5);
	    instance->clone()->scale(0.5);
		instance->scale(1.0);
	  });
	//.def("idot", &CustomLA::idot)
	//.def("copy_data_from", &CustomLA::copy_data_from);

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
         [](ROL::Algorithm<double> &instance, const std::string& str, const std::map<std::string, std::string>& params)
         {
           Teuchos::ParameterList parlist;
           parlist.sublist("Step").sublist("Line Search")
             .sublist("Descent Method").set("Type", "Quasi-Newton Method");
           parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
           parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
           parlist.sublist("Status Test").set("Iteration Limit",100);

           new (&instance) ROL::Algorithm<double>(str, parlist);
         })
    //.def("run", (std::vector<std::string> (ROL::Algorithm<double>::*)(ROL::Vector<double>&, ROL::Objective<double>&, bool, std::ostream&)) &ROL::Algorithm<double>::run);
    .def("run", [](ROL::Algorithm<double> &instance, ROL::Vector<double>& x, ROL::Objective<double>& obj)
	  {
	    instance.run(x, obj, true, std::cout);
	  });

  return m.ptr();
}
