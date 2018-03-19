#include <pybind11/stl.h>
namespace py = pybind11;

#include <ROL_ParameterList.hpp>

void dictToParameterList(py::dict param_dict, ROL::ParameterList& parlist) {
    for (auto item : param_dict) {
        auto key = item.first;
        auto value = item.second;

        if (py::isinstance<py::str>(key)) {
            auto key_ = std::string(py::str(key));
            if (py::isinstance<py::str>(value))
                parlist.set(key_, std::string(py::str(value)));
            else if (py::isinstance<py::bool_>(value))
                parlist.set(key_, value.cast<bool>());
            else if (py::isinstance<py::int_>(value))
                parlist.set(key_, value.cast<int>());
            else if (py::isinstance<py::float_>(value))
                parlist.set(key_, value.cast<double>());
            else if (py::isinstance<py::dict>(value)) {
                auto& sublist = parlist.sublist(std::string(py::str(key)));
                dictToParameterList(value.cast<py::dict>(), sublist);
            } else
                throw std::runtime_error(
                    std::string("Don't know what to do with value."));
        } else
            throw std::runtime_error(
                std::string("Don't know what to do with key."));
    }
}

void init_parameterlist(py::module& m) {
  py::class_<ROL::ParameterList, std::shared_ptr<ROL::ParameterList>>(
      m, "ParameterList", "Create a ParameterList object from an XML string")
      .def(py::init([](std::string xml_file) {
        return new ROL::ParameterList(
            *(ROL::getParametersFromXmlFile(xml_file)));
      }))
      .def(py::init([](py::dict param_dict, std::string name) {
        auto res = new ROL::ParameterList(name);
        dictToParameterList(param_dict, *res);
        return res;
      }))
      .def("print", [](ROL::ParameterList& instance) { instance.print(); });
}
