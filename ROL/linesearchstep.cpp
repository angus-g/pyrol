#include <pybind11/pybind11.h>
namespace py = pybind11;
#include "py_shared_ptr.hpp"
PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

#include <ROL_LineSearchStep.hpp>
#include <ROL_Secant.hpp>

class InitBFGS : public ROL::Secant<double> {
public:
  using ROL::Secant<double>::Secant;

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

  }
};

class PyInitBFGS : public InitBFGS {
public:
  using InitBFGS::InitBFGS;

  void applyH0(ROL::Vector<double> &Hv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE(void, ROL::Secant<double>, applyH0, Hv, v);
  }
};

void init_linesearchstep(py::module& m) {
  py::class_<ROL::Secant<double>, std::shared_ptr<ROL::Secant<double>>>(m, "Secant");

  // class, trampoline, reference type and parent class
  py::class_<InitBFGS, PyInitBFGS, py_shared_ptr<InitBFGS>, ROL::Secant<double>>(m, "InitBFGS")
    .def(py::init<>())
    .def("applyH0", &ROL::Secant<double>::applyH0);

  py::class_<ROL::Step<double>, std::shared_ptr<ROL::Step<double>>>(m, "Step");

  py::class_<ROL::LineSearchStep<double>, std::shared_ptr<ROL::LineSearchStep<double>>, ROL::Step<double>>(m, "LineSearchStep")
    .def(py::init([](ROL::ParameterList &params, std::shared_ptr<ROL::Secant<double>> &secant) {
      return new ROL::LineSearchStep<double>(params, nullptr, secant);
    }));
}
