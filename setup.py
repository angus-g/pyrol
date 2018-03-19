import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

class TrilinosInstallation(object):

    def __init__(self, user=False):
        self.user = user

    def get_path(self):
        import roltrilinos
        return roltrilinos.get_trilinos_path()

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

def set_rpath_linux(library, rpath):
    subprocess.check_call(["patchelf", "--set-rpath", rpath, library])

# def set_rpath_mac(library, rpath):
#     subprocess.check_call(["install_name_tool", "-add_rpath", rpath, library])

def fix_rpaths(trilinos_install_path):
    lib_path=os.path.join(trilinos_install_path, "lib")
    try:
        for s in os.listdir(lib_path):
            if ".so" in s:
                set_rpath_linux(os.path.join(lib_path, s), lib_path)
            # if ".dylib" in s
            #     set_rpath_mac(os.path.join(lib_path, s), lib_path)
    except:
        print("Could not set rpath; this might lead to a broken installation.")
        pass


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir + '\\' + ext.name,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]


        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']
        trilinos_installation = TrilinosInstallation()
        trilinos_dir = trilinos_installation.get_path()
        try:
            fix_rpaths(trilinos_dir)
        except:
            pass

        cmake_args += ['-DTRILINOS_DIR:Path=' + trilinos_dir]
        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='ROL',
    version='0.0.9',
    author='Chris Richardson, Greg von Winckel, Florian Wechsung',
    author_email='wechsung@maths.ox.ac.uk',
    description='A python wrapper for the ROL package.',
    long_description='',
    ext_modules=[CMakeExtension('ROL')],
    packages=['ROL'],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
