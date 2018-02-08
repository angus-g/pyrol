#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <ROL_Algorithm.hpp>
#include <ROL_Bounds.hpp>

void init_algorithm(py::module& m) {
    //
    // ROL::Algorithm<double>
    //
    py::class_<ROL::Algorithm<double>>(m, "Algorithm")
        .def(py::init<const std::string&, ROL::ParameterList&>())
        .def("run",
             [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
                ROL::Objective<double>& obj) {
                 instance.run(x, obj, true, std::cout);
             })
        .def("run",
             [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
                ROL::Objective<double>& obj, ROL::Bounds<double>& bnd) {
                 instance.run(x, obj, bnd, true, std::cout);
             })
        .def("run",
             [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
                ROL::Vector<double>& l, ROL::Objective<double>& obj,
                ROL::Constraint<double>& con) {
                 instance.run(x, l, obj, con, true, std::cout);
             })
        .def("run",
             [](ROL::Algorithm<double>& instance, ROL::Vector<double>& x,
                ROL::Vector<double>& l, ROL::Objective<double>& obj,
                ROL::Constraint<double>& con, ROL::Bounds<double>& bnd) {
                 instance.run(x, l, obj, con, bnd, true, std::cout);
             })
        .def("run",
             [](ROL::Algorithm<double>& instance,
                ROL::OptimizationProblem<double>& opt) {
                 instance.run(opt, true, std::cout);
             })
        .def("get_state", [](ROL::Algorithm<double>& instance) {
            return instance.getState();
        });
}
