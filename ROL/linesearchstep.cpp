#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_LineSearchStep.hpp>
#include <ROL_Secant.hpp>

class PySecant : public ROL::Secant<double> {
public:
  using ROL::Secant<double>::Secant;

  void applyH(ROL::Vector<double> &Hv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE_PURE(void, ROL::Secant<double>, applyH, Hv, v);
  }

  void applyB(ROL::Vector<double> &Bv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE_PURE(void, ROL::Secant<double>, applyB, Bv, v);
  }

  void applyH0(ROL::Vector<double> &Hv, const ROL::Vector<double> &v) const override {
    PYBIND11_OVERRIDE(void, ROL::Secant<double>, applyH0, Hv, v);
  }
};

class PyLineSearch : public ROL::LineSearch<double> {
public:
  using ROL::LineSearch<double>::LineSearch;

  void run(double &alpha, double &fval, int &ls_neval, int &ls_ngrad,
	   const double &gs, const ROL::Vector<double> &s, const ROL::Vector<double> &x,
	   ROL::Objective<double> &obj, ROL::BoundConstraint<double> &con) override {
    PYBIND11_OVERRIDE_PURE(void, ROL::LineSearch<double>, run, alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con);
  }
};

void init_linesearchstep(py::module& m) {
  py::class_<ROL::Secant<double>, PySecant>(m, "Secant")
    .def(py::init<>())
    .def("applyH0", &ROL::Secant<double>::applyH0);

  py::class_<ROL::LineSearch<double>, PyLineSearch>(m, "LineSearch")
    .def(py::init<ROL::ParameterList&>());

  py::class_<ROL::LineSearchStep<double>>(m, "LineSearchStep")
    .def(py::init<ROL::ParameterList&,
	 std::shared_ptr<ROL::LineSearch<double>>,
	 std::shared_ptr<ROL::Secant<double>>>(),
	 py::arg("parlist"),
	 py::arg("linesearch") = (ROL::LineSearch<double>*)nullptr,
	 py::arg("secant") = (ROL::Secant<double>*)nullptr);
}
