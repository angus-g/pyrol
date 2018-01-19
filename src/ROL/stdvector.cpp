#include <pybind11/stl.h>
namespace py = pybind11;

#include <ROL_StdVector.hpp>

void init_stdvector(py::module& m) {
  //
  // ROL::StdVector<double>
  //
  py::class_<ROL::StdVector<double>, std::shared_ptr<ROL::StdVector<double>>,
             ROL::Vector<double>>(m, "StdVector")
      .def(py::init([](int n) {
        auto tp = std::make_shared<std::vector<double>>(n, 0.0);
        auto ptr = std::make_shared<ROL::StdVector<double>>(tp);
        return ptr;
      }))
      .def("norm", &ROL::StdVector<double>::norm, "L2 norm of the vector")
      .def("dimension", &ROL::StdVector<double>::dimension,
           "Size of the vector")
      .def("__setitem__",
           [](ROL::StdVector<double>& vec, const int& idx, const double& val) {
             auto vvec = vec.getVector();
             if (idx >= (int)vvec->size())
               throw py::index_error();
             else
               (*vvec)[idx] = val;
           })
      .def("__getitem__",
           [](ROL::StdVector<double>& vec, const py::slice& slice) {
             auto vvec = vec.getVector();
             py::size_t start, stop, step, slicelength;
             if (!slice.compute((py::size_t)(vvec->size()), &start, &stop,
                                &step, &slicelength))
               throw py::error_already_set();
             auto res = std::make_shared<std::vector<double>>();
             res->resize(slicelength);
             for (unsigned int i = start; i < stop; i = i + step) {
               res->push_back((*vvec)[i]);
             }
           })
      .def("__getitem__",
           [](ROL::StdVector<double>& vec, const int& idx) {
             auto vvec = vec.getVector();
             if (idx >= (int)vvec->size()) {
               throw py::index_error();
             } else {
               return (*vvec)[idx];
             }
           })
      .def("scale", &ROL::StdVector<double>::scale,
           "Multiply the vector by a scalar")
      .def("checkVector",
           [](std::shared_ptr<ROL::StdVector<double>>& x,
              std::shared_ptr<ROL::StdVector<double>>& y,
              std::shared_ptr<ROL::StdVector<double>>& z)
               -> std::vector<double> {
             return x->checkVector(*y, *z, true, std::cout);
           },
           "Check the accuracy of the linear algebra implementation");
}
