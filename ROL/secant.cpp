#include <pybind11/pybind11.h>
namespace py = pybind11;
#include "py_shared_ptr.hpp"
PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

#include "vector.hpp"
#include <ROL_Types.hpp>

#include <ROL_Secant.hpp>

#include <fstream>
#define FMT_HEADER_ONLY
#include <fmt/core.h>

class InitBFGS : public ROL::Secant<double> {
public:
  using ROL::Secant<double>::Secant;

  virtual void scaleH0(ROL::Vector<double> &Hv) const {
    const auto& state = ROL::Secant<double>::get_state();

    if (state->iter != 0 && state->current != -1) {
      double yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
      Hv.scale(state->product[state->current] / yy);
    }
  }

  virtual void scaleB0(ROL::Vector<double> &Bv) const {
    const auto& state = ROL::Secant<double>::get_state();

    if (state->iter != 0 && state->current != -1) {
      double yy = state->gradDiff[state->current]->dot(*(state->gradDiff[state->current]));
      Bv.scale(yy / state->product[state->current]);
    }
  }

  void applyH(ROL::Vector<double> &Hv, const ROL::Vector<double> &v) const override {
    const auto& state = ROL::Secant<double>::get_state();
    double zero(0);

    Hv.set(v.dual());
    std::vector<double> alpha(state->current+1, zero);
    for (int i = state->current; i >= 0; i--) {
      alpha[i] = state->iterDiff[i]->dot(Hv);
      alpha[i] /= state->product[i];

      Hv.axpy(-alpha[i], (state->gradDiff[i])->dual());
    }

    auto tmp = Hv.clone();
    applyH0(*tmp, Hv.dual());
    Hv.set(*tmp);

    double beta(0);
    for (int i = 0; i <= state->current; i++) {
      beta = Hv.dot((state->gradDiff[i])->dual());
      beta /= state->product[i];

      Hv.axpy(alpha[i] - beta, *(state->iterDiff[i]));
    }
  }

  void applyB(ROL::Vector<double> &Bv, const ROL::Vector<double> &v) const override {
    const auto& state = ROL::Secant<double>::get_state();
    double one(1);

    applyB0(Bv, v);

    std::vector<std::shared_ptr<ROL::Vector<double>>> a(state->current + 1), b(state->current + 1);
    double bv(0), av(0), bs(0), as(0);
    for (int i = 0; i <= state->current; i++) {
      b[i] = Bv.clone();
      b[i]->set(*(state->gradDiff[i]));
      b[i]->scale(one / sqrt(state->product[i]));
      bv = v.dot(b[i]->dual());
      Bv.axpy(bv, *b[i]);

      a[i] = Bv.clone();
      applyB0(*a[i], *(state->iterDiff[i]));

      for (int j = 0; j < i; j++) {
	bs = (state->iterDiff[i])->dot(b[j]->dual());
	a[i]->axpy(bs, *b[j]);

	as = (state->iterDiff[i])->dot(a[j]->dual());
	a[i]->axpy(-as, *a[j]);
      }

      as = (state->iterDiff[i])->dot(a[i]->dual());
      a[i]->scale(one / sqrt(as));

      av = v.dot(a[i]->dual());
      Bv.axpy(-av, *a[i]);
    }
  }
};

class PyInitBFGS : public InitBFGS {
public:
  using InitBFGS::InitBFGS;

  void scaleH0(ROL::Vector<double> &Hv) const override {
    PYBIND11_OVERRIDE(void, InitBFGS, scaleH0, Hv);
  }

  void scaleB0(ROL::Vector<double> &Bv) const override {
    PYBIND11_OVERRIDE(void, InitBFGS, scaleB0, Bv);
  }

  void applyH0(ROL::Vector<double> &Hv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE(void, InitBFGS, applyH0, Hv, v);
  }

  void applyB0(ROL::Vector<double> &Bv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE(void, InitBFGS, applyB0, Bv, v);
  }
};

void serialise_secant(const ROL::Secant<double> &v, int rank) {
  std::ofstream os(fmt::format("rol_secant_{}.cereal", rank), std::ios::binary);
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(v);
}

void load_secant(ROL::Secant<double> &v, int rank) {
  std::ifstream is(fmt::format("rol_secant_{}.cereal", rank), std::ios::binary);
  cereal::BinaryInputArchive iarchive(is);

  iarchive(v);
}

void init_secant(py::module& m) {
  m.def("serialise_secant", &serialise_secant, "Serialise a ROL Secant using Cereal",
	py::arg("secant"), py::arg("rank") = 0);
  m.def("load_secant", &load_secant, "Load a ROL Secant from Cereal archive",
	py::arg("secant"), py::arg("rank") = 0);

  py::class_<ROL::Secant<double>, std::shared_ptr<ROL::Secant<double>>>(m, "Secant");

  // class, trampoline, reference type and parent class
  py::class_<InitBFGS, PyInitBFGS, py_shared_ptr<InitBFGS>, ROL::Secant<double>>(m, "InitBFGS")
    .def(py::init<int>())
    .def("applyH0", &InitBFGS::applyH0)
    .def("applyB0", &InitBFGS::applyB0)
    .def("scaleH0", &InitBFGS::scaleH0)
    .def("scaleB0", &InitBFGS::scaleB0);
}
