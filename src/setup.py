import os
from setuptools import setup, Extension

tri_dir=os.environ["TRILINOS_DIR"]
if not os.path.exists(tri_dir+"/include"):
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
                    include_dirs=["pybind11/include",
                                  tri_dir+"/include",
                                  "/usr/include/mpi"],
                    library_dirs=[tri_dir+"/lib"],
                    libraries=["teuchoscore", "teuchosnumerics", "teuchosparameterlist", "teuchoscomm"],
                    extra_compile_args=["-std=c++11"]
                    )],
      )
