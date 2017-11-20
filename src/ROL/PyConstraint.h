#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Constraint.hpp>

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

    virtual void update(const ROL::Vector<double> &x, bool flag = true,
                        int iter = -1) override {
        PYBIND11_OVERLOAD(void, ROL::Constraint<double>, update, x, flag, iter);
    }
};
