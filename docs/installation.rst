Installation
============

We assume that you have Trilinos as well as eigen installed on your machine.

.. code-block:: bash

    git clone https://bitbucket.org/pyrol/pyrol/
    cd src
    git clone https://github.com/pybind/pybind11

We use cmake for the configuration. Sometime there are issues with finding the correct python version.
On macOS, try running 

.. code-block:: bash

    cmake -DPYTHON_LIBRARY=$(python-config --prefix)/lib/libpython2.7.dylib \-DPYTHON_INCLUDE_DIR=$(python-config --prefix)/include/python2.7 -DTRILINOS_DIR="/path/to/Trilinos/" ..

on Ubuntu run

.. code-block:: bash

    cmake -DTRILINOS_DIR="/path/to/Trilinos/" -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python2.7 -DPYTHON_INCLUDE_DIR:PATH=/usr/include/python2.7 -DPYTHON_LIBRARY:FILEPATH=/usr/lib/x86_64-linux-gnu/libpython2.7.so ..


