#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Constraint.hpp>
#include "py_shared_ptr.hpp"
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, py_shared_ptr<T>);

class PyConstraint : public ROL::Constraint<double> {
   public:
    void value(ROL::Vector<double> &c, const ROL::Vector<double> &x,
               double &tol) override {
        PYBIND11_OVERLOAD_PURE(void, ROL::Constraint<double>, value, c, x, tol);
    }

    virtual void applyJacobian(ROL::Vector<double> &jv,
                               const ROL::Vector<double> &v,
                               const ROL::Vector<double> &x,
                               double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Constraint<double>, applyJacobian, jv, v,
                          x, tol);
    }

    virtual void applyAdjointJacobian(ROL::Vector<double> &ajv,
                                      const ROL::Vector<double> &v,
                                      const ROL::Vector<double> &x,
                                      double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Constraint<double>, applyAdjointJacobian,
                          ajv, v, x, tol);
    }
    virtual void applyAdjointHessian(ROL::Vector<double> &ahuv,
                                     const ROL::Vector<double> &u,
                                     const ROL::Vector<double> &v,
                                     const ROL::Vector<double> &x, double &tol) override {
        PYBIND11_OVERLOAD(void, ROL::Constraint<double>, applyAdjointHessian,
            ahuv, u, v, x, tol);
      }


    virtual void update(const ROL::Vector<double> &x, bool flag = true,
                        int iter = -1) override {
        PYBIND11_OVERLOAD(void, ROL::Constraint<double>, update, x, flag, iter);
    }
};

void init_constraint(py::module &m) {
    //
    // ROL::Constraint<double>
    //
    py::class_<ROL::Constraint<double>, PyConstraint,
               py_shared_ptr<ROL::Constraint<double>>>(m, "Constraint")
        .def(py::init<>())
        //.def("value", &ROL::Constraint<double>::value)
        //.def("applyJacobian", &ROL::Constraint<double>::applyJacobian)
        //.def("applyAdjointJacobian", [](ROL::Constraint<double>& instance,
        //  ROL::Vector<double> &ajv, const ROL::Vector<double> &v,
        //  const ROL::Vector<double> &x, double &tol)
        //  {
        //    instance.applyAdjointJacobian(ajv, v, x, tol);
        //  })
        .def("checkApplyJacobian",
             [](ROL::Constraint<double> &instance, ROL::Vector<double> &x,
                ROL::Vector<double> &v, ROL::Vector<double> &jv, int steps,
                int order) {
                 return instance.checkApplyJacobian(x, v, jv, true, std::cout, steps,
                                             order);
             })
        .def("checkAdjointConsistencyJacobian",
             [](ROL::Constraint<double> &instance, ROL::Vector<double> &w,
                ROL::Vector<double> &v, ROL::Vector<double> &x) {
                 return instance.checkAdjointConsistencyJacobian(w, v, x, true,
                                                          std::cout);
             })
        .def("checkApplyAdjointHessian",
             [](ROL::Constraint<double> &instance, ROL::Vector<double> &x,
                ROL::Vector<double> &u, ROL::Vector<double> &v, 
                 ROL::Vector<double> &hv, int numSteps, int order){
                 return instance.checkApplyAdjointHessian(x, u, v, hv, true, std::cout, numSteps, order);
             });
}
