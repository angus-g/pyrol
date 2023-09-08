#include <pybind11/pybind11.h>
namespace py = pybind11;


void init_algorithmstate(py::module&);
void init_bounds(py::module&);
void init_constraint(py::module&);
void init_objective(py::module&);
void init_stdvector(py::module&);
void init_vector(py::module&);
void init_statustest(py::module&);
void init_optimizationproblem(py::module&);
void init_optimizationsolver(py::module&);
void init_secant(py::module&);
void init_algorithm(py::module&);
void init_problem(py::module&);
void init_parameterlist(py::module&);

PYBIND11_MODULE(_ROL, m) {
  m.doc() =
      "PyROL provides Python wrappers for a subset of the"
      "Trilinos ROL library.";
  m.attr("__version__") = "0.2.dev0";

  init_vector(m);
  init_algorithmstate(m);
  init_bounds(m);
  init_constraint(m);
  init_objective(m);
  init_stdvector(m);
  init_statustest(m);
  init_optimizationproblem(m);
  init_optimizationsolver(m);
  init_secant(m);
  init_algorithm(m);
  init_problem(m);
  init_parameterlist(m);
}
