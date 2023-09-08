#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>

#include <ROL_Objective.hpp>

class PyObjective : public ROL::Objective<double> {
   public:
    double value(const ROL::Vector<double> &x, double &tol) override {
        PYBIND11_OVERLOAD_PURE(double, ROL::Objective<double>, value, x, tol);
    }
    void gradient(ROL::Vector<double> &g, const ROL::Vector<double> &x,
                  double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Objective<double>, gradient, g, x, tol);
    }

    void hessVec(ROL::Vector<double> &hv, const ROL::Vector<double> &v,
                 const ROL::Vector<double> &x, double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Objective<double>, hessVec, hv, v, x, tol);
    }

    void invHessVec(ROL::Vector<double> &hv, const ROL::Vector<double> &v,
                    const ROL::Vector<double> &x, double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Objective<double>, invHessVec, hv, v, x,
                          tol);
    }

    virtual void update(const ROL::Vector<double> &x, bool flag = true,
                        int iter = -1) override {
        PYBIND11_OVERLOAD(void, ROL::Objective<double>, update, x, flag, iter);
    }

  virtual void update(const ROL::Vector<double> &x, ROL::UpdateType type, int iter = -1) override {
    PYBIND11_OVERRIDE(void, ROL::Objective<double>, update, x, type, iter);
  }
};

void init_objective(py::module &m) {
  //
  // ROL::Objective<double>
  //
  py::class_<ROL::Objective<double>, PyObjective, std::shared_ptr<ROL::Objective<double>>>
      objective(
          m, "Objective",
          "Base class for the objective class. Python objectives need to"
          "inherit from this class.");
  objective.def(py::init<>())
  .def("checkGradient", [](ROL::Objective<double> &instance, ROL::Vector<double>& x,
       std::shared_ptr<ROL::Vector<double>> d, int steps, int order, double scale)
      {
        if(d == nullptr)
        {
          d = x.dual().clone();
          double tol = 0.0;
          instance.gradient(*d, x, tol);
          d->scale(scale);
          return instance.checkGradient(x, d->dual(), true, std::cout, steps, order);
        }
        else {
          d->scale(scale);
          return instance.checkGradient(x, *d, true, std::cout, steps, order);
        }
      }, py::arg("x"), py::arg("d")=(ROL::Vector<double>*)nullptr, py::arg("steps")=4,
      py::arg("order")=1, py::arg("scale")=1.0
      )
    .def("checkHessVec",
        [](ROL::Objective<double> &instance, ROL::Vector<double> &x,
          std::shared_ptr<ROL::Vector<double>> v, int steps, int order, double scale) {
        if(v == nullptr)
        {
          v = x.dual().clone();
          double tol = 0.0;
          instance.gradient(*v, x, tol);
          v->scale(scale);
          return instance.checkHessVec(x, v->dual(), true, std::cout, steps, order);
        }
        else {
          v->scale(scale);
          return instance.checkHessVec(x, *v, true, std::cout, steps, order);
        }
        }, py::arg("x"), py::arg("v")=(ROL::Vector<double>*)nullptr, py::arg("steps")=4,
        py::arg("order")=1, py::arg("scale")=1.0);

  py::enum_<ROL::UpdateType>(m, "UpdateType")
    .value("Initial", ROL::UpdateType::Initial)
    .value("Accept", ROL::UpdateType::Accept)
    .value("Revert", ROL::UpdateType::Revert)
    .value("Trial", ROL::UpdateType::Trial)
    .value("Temp", ROL::UpdateType::Temp);
}

CEREAL_REGISTER_TYPE(PyObjective);
CEREAL_REGISTER_POLYMORPHIC_RELATION(ROL::Objective<double>, PyObjective);
