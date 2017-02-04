#include <pybind11/pybind11.h>
#include <Python.h>
namespace py = pybind11;

#include <ROL_Vector.hpp>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"
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
				// A wiser man than me maybe knows why this ref count increase is necessary...
				Py_INCREF(res.ptr());
				return Teuchos::rcp(res.cast<std::shared_ptr<CustomLA>>());
			}
			else
				py::pybind11_fail("Tried to call pure virtual function 'clone'.");
		}

};
