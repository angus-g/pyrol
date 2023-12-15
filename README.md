pyrol
=====

pyrol is a wrapper for the Trilinos
[ROL](https://trilinos.github.io/rol) library using
[pybind11](https://github.com/pybind/pybind11). It also supports
serialisation of in-progress optimisation using [cereal](https://github.com/USCiLab/cereal).

Installation
------------

pyrol will build a compatible version of the Trilinos shared library
as part of its build process. For this, a C++ compiler, CMake, BLAS
and LAPACK implementations are all required. With these components in
place, installation should be as simple as `pip install`. For more
complicated environments, CMake arguments can be provided through the
`CMAKE_ARGS` environment variable.

License
-------

PyROL is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.  PyROL is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU Lesser General Public License for more details.  You should
have received a copy of the GNU Lesser General Public License along
with PyROL. If not, see <http://www.gnu.org/licenses/>.
