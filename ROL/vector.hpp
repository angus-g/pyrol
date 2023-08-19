// register polymorphic serialisation to binary archive
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/memory.hpp>

#include <ROL_Vector.hpp>

#include <iostream>

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
        PYBIND11_OVERLOAD(int, ROL::Vector<double>, dimension, );
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
      std::function<double(double, double)> stdr = [&r](double a, double r0) {
        r.reduce(a, r0);
        return r0;
      };
      PYBIND11_OVERLOAD_PURE(double, ROL::Vector<double>, reduce,
                             py::cpp_function(stdr), r.initialValue());
    }

    void applyUnary(const ROL::Elementwise::UnaryFunction<double>& f) override {
        std::function<double(double)> stdf = [&f](double a){return f.apply(a);};
        PYBIND11_OVERLOAD_PURE(void, ROL::Vector<double>, applyUnary, py::cpp_function(stdf));
    }

    void applyBinary(const ROL::Elementwise::BinaryFunction<double>& f,
                     const ROL::Vector<double>& x) override {
        if (dimension() != x.dimension()) {
            throw std::length_error(
                "Error: Vectors must have the same dimension.");
        }
        std::function<double(double, double)> stdf = [&f](double a, double b){return f.apply(a, b);};
        PYBIND11_OVERLOAD_PURE(void, ROL::Vector<double>, applyBinary, py::cpp_function(stdf), x);
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


template <class Archive>
void save(Archive &archive, PyVector const &v) {
  py::gil_scoped_acquire gil;

  auto pickle = py::module_::import("pickle");
  auto pickled = pickle.attr("dumps")(py::cast(v));
  auto result = pickled.cast<std::string>();

  archive(result);
}

// we need this so that the polymorphism business knows how to deal with PyVector
template <class Archive>
void load(Archive &archive, PyVector &v) {
}

template <class Archive>
void load(Archive &archive, cereal::memory_detail::PtrWrapper<std::shared_ptr<PyVector> &> &wrapper) {
  uint32_t id;
  archive(CEREAL_NVP_("id", id));

  if (id & cereal::detail::msb_32bit) {
    py::gil_scoped_acquire gil;

    auto pickle = py::module_::import("pickle");
    std::string pickled;
    archive(CEREAL_NVP_("data", pickled));
    auto loaded = pickle.attr("loads")(py::bytes(pickled));

    auto keep_alive = std::make_shared<py::object>(loaded);
    auto p = loaded.cast<PyVector *>();
    std::shared_ptr<PyVector> ptr(keep_alive, p);

    archive.registerSharedPointer(id, ptr);
    wrapper.ptr = std::move(ptr);
  } else {
    wrapper.ptr = std::static_pointer_cast<PyVector>(archive.getSharedPointer(id));
  }
}

CEREAL_REGISTER_TYPE(PyVector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(ROL::Vector<double>, PyVector);
