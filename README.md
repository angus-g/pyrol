README
------

We are wrapping [ROL](trilinos.org/packages/rol/) with pybind11. Please see [Issues](../../issues) to contribute proposals, tasks, bugs etc.


Florian Wechsung (@florianwechsung)

Chris Richardson (@chris_richardson)

INSTALLATION
------------

PyROL needs to link to an installation of ROL which is part of the [Trilinos](https://github.com/trilinos/Trilinos) package.
Since this is a C++ package, we need to do some rpath magic to make things work. If you're on linux you should install `patchelf` first, e.g. by 

    sudo apt install patchelf

Running 
    
    pip3 install roltrilinos

will then install Trilinos with the proper options set.
You can then either

    pip3 install ROL

or if you want to contribute to ROL, you simply run

    git clone git@bitbucket.org:pyrol/pyrol.git
    cd pyrol
    git submodule update --init
    pip3 install -e .


License
-------

PyROL is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PyROL is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with PyROL. If not, see <http://www.gnu.org/licenses/>.
