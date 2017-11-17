#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Constraint.hpp>

class PyConstraint : public ROL::Constraint<double>
{
	public:

		void value(ROL::Vector<double>& c, const ROL::Vector<double>& x, double &tol ) override
		{
			py::gil_scoped_acquire gil;
			py::function overload = py::get_overload(this, "value");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(c, x, tol).cast<void>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'value'.");
		}

		virtual void applyJacobian(ROL::Vector<double> &jv, const ROL::Vector<double> &v,
				const ROL::Vector<double> &x, double &tol) override
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "applyJacobian");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(jv, v, x, tol).cast<void>();
			else
				ROL::Constraint<double>::applyJacobian(jv, v, x, tol);
		}

		virtual void applyAdjointJacobian(ROL::Vector<double> &ajv, const ROL::Vector<double> &v,
				const ROL::Vector<double> &x, double &tol) override
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "applyAdjointJacobian");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(ajv, v, x, tol).cast<void>();
			else
				ROL::Constraint<double>::applyAdjointJacobian(ajv, v, x, tol);
		}

		virtual void update( const ROL::Vector<double> &x, bool flag = true, int iter = -1 ) override
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "update");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(x, flag, iter).cast<void>();
			else
				ROL::Constraint<double>::update(x, flag, iter);
		}
};
