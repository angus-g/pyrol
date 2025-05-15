cmake_minimum_required(VERSION 3.14)
project(_ROL VERSION 0.4.0)

# Required for std::filesystem
set(CMAKE_CXX_STANDARD 17)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG)

pybind11_add_module(_ROL MODULE ROL/ROL.cpp
  ROL/algorithmstate.cpp
  ROL/augmentedlagrangian.cpp
  ROL/bounds.cpp
  ROL/constraint.cpp
  ROL/objective.cpp
  ROL/statustest.cpp
  ROL/optimizationproblem.cpp
  ROL/optimizationsolver.cpp
  ROL/secant.cpp
  ROL/algorithm.cpp
  ROL/problem.cpp
  ROL/parameterlist.cpp
  ROL/stdvector.cpp
  ROL/vector.cpp)

#target_include_directories(_ROL PRIVATE ${ROL_INCLUDE_DIRS})
target_link_libraries(_ROL PRIVATE rol)
# set_target_properties(_ROL PROPERTIES COMPILE_DEFINITIONS "ENABLE_PYROL")
