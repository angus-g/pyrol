import os, subprocess
from setuptools import setup, Extension

# Detect pybind11 installation and add include path
from setuptools.command.build_ext import build_ext
class build_ext_pybind11(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
        import pybind11
        self.include_dirs.append(pybind11.get_include())

# Try to get the MPI flags
try:
    mpicxx = subprocess.check_output(["mpicxx","-show"])
    ex_flags = mpicxx.strip()
    ex_flags = ex_flags.split(" ")[1:]
except IOError:
    print 'Cannot query mpicxx compiler used to compile Trilinos\n'
    print 'Make sure that mpicxx exists, and returns flags with "mpicxx -show"'
    quit()

# Make sure we are using c++11
ex_flags.append("-std=c++11")

if "TRILINOS_DIR" not in os.environ:
    print "Please set the TRILINOS_DIR environment variable to point to the Trilinos"
    print "installation directory containing ROL and Teuchos\n"
    quit()

trilinos_dir = os.environ["TRILINOS_DIR"]
if not os.path.exists(trilinos_dir+"/include/ROL_config.h"):
    print "Cannot find Trilinos include directory with ROL. Make sure TRILINOS_DIR is set correctly.\n"
    quit()

setup(name="PyROL",
      version="0.1.1",
      author="Chris Richardson and Florian Wechsung",
      author_email="pyrol-dev@googlemail.com",
      license="LGPLv3",
      packages=['ROLUtils'],
      install_requires=["pybind11"],
      cmdclass={'build_ext': build_ext_pybind11},
      ext_modules=[
          Extension("ROL",
                    ["ROL/ROL.cpp"],
                    include_dirs = [trilinos_dir+"/include"],
                    library_dirs = [trilinos_dir+"/lib"],
                    runtime_library_dirs = [trilinos_dir+"/lib"],
                    libraries = ["teuchoscore", "teuchosnumerics", "teuchosparameterlist", "teuchoscomm"],
                    extra_compile_args = ex_flags,
                 )],
      )
