from __future__ import print_function
import os, sys, subprocess
from setuptools import setup, Extension

if "TRILINOS_DIR" not in os.environ:
    print("Please set the TRILINOS_DIR environment variable to point to the Trilinos")
    print("installation directory containing ROL and Teuchos\n")
    sys.exit(1)

trilinos_dir = os.environ["TRILINOS_DIR"]

if (sys.platform == 'darwin'):
    os.environ.setdefault('LDFLAGS','')
    os.environ['LDFLAGS'] += os.path.expandvars(' -Wl,-rpath,$TRILINOS_DIR/lib')

try:
    import pybind11
except ImportError:
    raise ImportError("Please pip install pybind11 first")
else:
    pybind11_get_include = pybind11.get_include()

# Try to get the MPI flags
try:
    mpicxx = subprocess.check_output(["mpicxx", "-show"])
    ex_flags = mpicxx.strip()
    ex_flags = ex_flags.decode().split(" ")[1:]
except IOError:
    print('Cannot query mpicxx compiler used to compile Trilinos\n')
    print('Make sure that mpicxx exists, and returns flags with "mpicxx -show"')
    sys.exit(1)

# Make sure we are using c++11
ex_flags.append("-std=c++11")

# Collect up libraries needed by mpi
ex_libs = [x[2:] for x in ex_flags if x[:2]=='-l']

if not os.path.exists(trilinos_dir+"/include/ROL_config.h"):
    print("Cannot find Trilinos include directory with ROL. Make sure TRILINOS_DIR is set correctly.\n")
    sys.exit(1)

ex_libs.extend(["teuchoscore", "teuchosnumerics", "teuchosparameterlist", "teuchoscomm"])

setup(name="PyROL",
      version="0.1.1",
      author="Chris Richardson and Florian Wechsung",
      author_email="pyrol-dev@googlemail.com",
      license="LGPLv3",
      packages=['ROL'],
      # setup_requires=["pybind11"],
      ext_modules=[
          Extension("_ROL",
                    ["ROL/ROL.cpp",
                     "ROL/algorithm.cpp",
                     "ROL/algorithmstate.cpp",
                     "ROL/bounds.cpp",
                     "ROL/constraint.cpp",
                     "ROL/moreauyosidapenalty.cpp",
                     "ROL/objective.cpp",
                     "ROL/optimizationproblem.cpp",
                     "ROL/optimizationsolver.cpp",
                     "ROL/parameterlist.cpp",
                     "ROL/stdvector.cpp",
                     "ROL/vector.cpp"],
                    define_macros = [("ROL_SHARED_POINTER", 1)],
                    include_dirs = [pybind11_get_include] + [trilinos_dir+"/include"],
                    library_dirs = [trilinos_dir+"/lib"],
                    runtime_library_dirs = [trilinos_dir+"/lib"],
                    libraries = ex_libs,
                    extra_compile_args = ex_flags,
                 )],
      )
