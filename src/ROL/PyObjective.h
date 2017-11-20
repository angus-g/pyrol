#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Objective.hpp>

class PyObjective : public ROL::Objective<double> {
   public:
    double value(const ROL::Vector<double> &x, double &tol) override {
        PYBIND11_OVERLOAD_PURE(double, ROL::Objective<double>, value, x, tol);
    }
    void gradient(ROL::Vector<double> &g, const ROL::Vector<double> &x,
                  double &tol) override {
        PYBIND11_OVERLOAD_PURE(void, ROL::Objective<double>, gradient, g, x,
                               tol);
    }
    virtual void update(const ROL::Vector<double> &x, bool flag = true,
                        int iter = -1) override {
        PYBIND11_OVERLOAD(void, ROL::Objective<double>, update, x, flag, iter);
    }
};
