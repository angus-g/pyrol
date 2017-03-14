#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_BoundConstraint.hpp>

class PyBoundConstraint : public ROL::BoundConstraint<double>
{
	public:

        PyBoundConstraint()
            : ROL::BoundConstraint<double>()
        {}

        PyBoundConstraint(const Teuchos::RCP<ROL::Vector<double> > &x_lo,
                const Teuchos::RCP<ROL::Vector<double> > &x_up,
                const double scale = 1)
            : ROL::BoundConstraint<double>(x_lo, x_up, scale)
        {}

		void project(ROL::Vector<double>& x) override
		{
			py::gil_scoped_acquire gil;
			py::function overload = py::get_overload(this, "project");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(x).cast<void>();
			else
				ROL::BoundConstraint<double>::project(x);
		}

};
