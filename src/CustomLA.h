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

		virtual int dimension() const {
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "dimension");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>().cast<int>();
			else
				py::pybind11_fail("Tried to call pure virtual function 'dimension'.");
		}

		Teuchos::RCP<ROL::Vector<double>> basis( const int i ) const {
			if(i >= this->dimension())
				throw py::index_error();
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "basis");
			if (overload)
			{
				auto res = overload.operator()<py::return_value_policy::reference>(i);
				Py_INCREF(res.ptr());
				return Teuchos::rcp(res.cast<std::shared_ptr<CustomLA>>());
			}
			else
				py::pybind11_fail("Tried to call pure virtual function 'basis'.");
		}


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

		virtual double reduce( const ROL::Elementwise::ReductionOp<double> &r ) const {
			double res = r.initialValue();
			for (int i = 0; i < this->dimension(); ++i) {
				r.reduce(this->getitem(i), res);
			}
			return res;
		}

        void applyUnary( const ROL::Elementwise::UnaryFunction<double> &f ) {
            uint dim  = dimension();
            for(uint i=0; i<dim; ++i) {
				setitem(i, f.apply(getitem(i)));
            }
        }

		void applyBinary( const ROL::Elementwise::BinaryFunction<double> &f, const ROL::Vector<double> &x ) {
			TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
					std::invalid_argument,
					"Error: Vectors must have the same dimension." );

			const CustomLA & ex = Teuchos::dyn_cast<const CustomLA>(x);
			for (uint i=0; i<dimension(); i++) {
				setitem(i, f.apply(getitem(i), ex.getitem(i)));
			}
		}

		virtual double getitem(const int& i) const
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "__getitem__");
			if (overload)
				return overload.operator()<py::return_value_policy::reference>(i).cast<double>();
			else
				py::pybind11_fail("Tried to call pure virtual function '__getitem__'.");
		}

		virtual void setitem(const int& i, const double& val) const
		{
			py::gil_scoped_acquire gil;
			pybind11::function overload = py::get_overload(this, "__setitem__");
			if (overload)
				overload.operator()<py::return_value_policy::reference>(i, val).cast<void>();
			else
				py::pybind11_fail("Tried to call pure virtual function '__setitem__'.");
		}


};
