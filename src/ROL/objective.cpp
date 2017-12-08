#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

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
};

void init_objective(py::module &m) {
    //
    // ROL::Objective<double>
    //
    py::class_<ROL::Objective<double>, PyObjective,
               std::shared_ptr<ROL::Objective<double>>>
        objective(
            m, "Objective",
            "Base class for the objective class. Python objectives need to"
            "inherit from this class.");
    objective.def(py::init<>())
        .def("checkGradient",
             [](ROL::Objective<double> &instance, ROL::Vector<double> &x,
                ROL::Vector<double> &d, int steps, int order) {
                 auto res = instance.checkGradient(x, d, true, std::cout, steps,
                                                   order);
                 return res;
             })
        .def("checkHessVec",
             [](ROL::Objective<double> &instance, ROL::Vector<double> &x,
                ROL::Vector<double> &v, int steps, int order) {
                 auto res =
                     instance.checkHessVec(x, v, true, std::cout, steps, order);
                 return res;
             });
}
