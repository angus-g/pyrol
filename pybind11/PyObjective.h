#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Objective.hpp>

class PyObjective : public ROL::Objective<double>
{
	public:

		double value( const ROL::Vector<double>& x, double &tol ) override
		{
			py::gil_scoped_acquire gil;
			py::function overload = py::get_overload(this, "value");
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
