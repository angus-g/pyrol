#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Vector.hpp>

class PyVector : public ROL::Vector<double> {
   public:
    virtual void plus(const ROL::Vector<double>& x) override {
        PYBIND11_OVERLOAD_PURE(void, ROL::Vector<double>, plus, x);
    }

    virtual void scale(const double alpha) override {
        PYBIND11_OVERLOAD_PURE(void, ROL::Vector<double>, scale, alpha);
    }

    virtual double dot(const ROL::Vector<double>& x) const override {
        PYBIND11_OVERLOAD_PURE(double, ROL::Vector<double>, dot, x);
    }

    virtual double norm() const override { return std::sqrt(this->dot(*this)); }

    virtual int dimension() const override {
        PYBIND11_OVERLOAD_PURE(int, ROL::Vector<double>, dimension, );
    }

    std::shared_ptr<ROL::Vector<double>> basis(const int i) const override {
        if (i >= this->dimension()) throw py::index_error();
        PYBIND11_OVERLOAD_PURE(std::shared_ptr<ROL::Vector<double>>,
                               ROL::Vector<double>, basis, i);
    }

    virtual std::shared_ptr<ROL::Vector<double>> clone() const override {
        // see https://github.com/pybind/pybind11/issues/1049
        auto self = py::cast(this);
        auto cloned = self.attr("clone")();

        auto keep_python_state_alive = std::make_shared<py::object>(cloned);
        auto ptr = cloned.cast<PyVector*>();

        // aliasing shared_ptr: points to `A_trampoline* ptr` but refcounts the
        // Python object
        return std::shared_ptr<ROL::Vector<double>>(keep_python_state_alive,
                                                    ptr);
    }

    virtual void axpy(const double alpha,
                      const ROL::Vector<double>& x) override {
        PYBIND11_OVERLOAD(void, ROL::Vector<double>, axpy, alpha, x);
    }

    virtual void zero() override {
        PYBIND11_OVERLOAD(void, ROL::Vector<double>, zero, );
    }

    virtual void set(const ROL::Vector<double>& x) override {
        PYBIND11_OVERLOAD(void, ROL::Vector<double>, set, x);
    }

    virtual double reduce(
        const ROL::Elementwise::ReductionOp<double>& r) const override {
        double res = r.initialValue();
        for (int i = 0; i < this->dimension(); ++i) {
            r.reduce(this->getitem(i), res);
        }
        return res;
    }

    void applyUnary(const ROL::Elementwise::UnaryFunction<double>& f) override {
        uint dim = dimension();
        for (uint i = 0; i < dim; ++i) {
            setitem(i, f.apply(getitem(i)));
        }
    }

    void applyBinary(const ROL::Elementwise::BinaryFunction<double>& f,
                     const ROL::Vector<double>& x) override {
        if (dimension() != x.dimension()) {
            throw std::length_error(
                "Error: Vectors must have the same dimension.");
        }
        const PyVector& ex = dynamic_cast<const PyVector&>(x);
        for (int i = 0; i < dimension(); ++i) {
            setitem(i, f.apply(getitem(i), ex.getitem(i)));
        }
    }

    virtual double getitem(const int& i) const {
        PYBIND11_OVERLOAD_PURE_NAME(double, ROL::Vector<double>, "__getitem__",
                                    getitem, i);
    }

    virtual void setitem(const int& i, const double& val) const {
        PYBIND11_OVERLOAD_PURE_NAME(void, ROL::Vector<double>, "__setitem__",
                                    setitem, i, val);
    }

    virtual void print(std::ostream& outStream) const override {
        PYBIND11_OVERLOAD(void, ROL::Vector<double>, print, outStream);
    }
};

void init_vector(py::module& m) {
    py::class_<ROL::Vector<double>, std::shared_ptr<ROL::Vector<double>>,
               PyVector>(m, "Vector")
        .def(py::init<>());
}
