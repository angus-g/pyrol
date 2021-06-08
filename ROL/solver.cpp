#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <cereal/archives/binary.hpp>

#include <map>

#include <ROL_StatusTest.hpp>
#include <ROL_Solver.hpp>

void init_solver(py::module &m) {
  py::class_<ROL::Solver<double>,
	     std::shared_ptr<ROL::Solver<double>>>(m, "Solver")
    .def(py::init<std::shared_ptr<ROL::Problem<double>>&, ROL::ParameterList&>())
    .def("solve",
	 [](ROL::Solver<double> &instance, std::shared_ptr<ROL::StatusTest<double>> &status, bool print_output) {
	   if (print_output)
	     instance.solve(std::cout, status);
	   else
	     instance.solve(status);
	 },
	 py::arg("status") = (ROL::StatusTest<double>*)nullptr,
	 py::arg("print_output") = true)
    .def("getAlgorithmState", &ROL::Solver<double>::getAlgorithmState);
}
