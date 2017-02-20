import os
from setuptools import setup, Extension

import pybind11

trilinos_dir=os.environ["TRILINOS_DIR"]
if not os.path.exists(trilinos_dir+"/include"):
    print "Cannot find Trilinos. Set TRILINOS_DIR"
    quit()

setup(name="PyROL",
      version="0.1.1",
      author="Chris Richardson and Florian Wechsung",
      author_email="pyrol-dev@googlemail.com",
      packages=["ROL"],
      license="LGPLv3",
      install_requires=[""],
      ext_modules=[
          Extension("ROL",
                    ["ROL/ROL.cpp"],
                    include_dirs = [pybind11.get_include(), trilinos_dir+"/include", "/usr/include/mpi"],
                    library_dirs = [trilinos_dir+"/lib"],
                    runtime_library_dirs = [trilinos_dir+"/lib"],
                    libraries = ["teuchoscore", "teuchosnumerics", "teuchosparameterlist", "teuchoscomm", "mpi_cxx"],
                    extra_compile_args = ["-std=c++11"]
                    )],
      )
