Introduction
============

ROL (Rapid Optimization Library) is an optimization library written in C++ and is part of Trilinos. PyROL is a Python wrapper for ROL using the <https://github.com/pybind/pybind11> library. It allows you to create objectives, constraints and vectors all in python. For the vectors you can use your own underlying data storage.

This enables the use of ROL for large scale PDE Constrained Optimization with FEniCS or Firedrake. 
Have a look at the examples to see how to to solve the motherproblem of PDE Constrained Optimization with either of the two libraries.
