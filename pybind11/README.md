# Installation on MAC

Installation on mac is tricky due to some CMake issues with mac-python vs homebrew-python.


	mkdir build; cd build
	cmake -DPYTHON_LIBRARY=$(python-config --prefix)/lib/libpython2.7.dylib -DPYTHON_INCLUDE_DIR=$(python-config --prefix)/include/python2.7 -DTRILINOS_DIR="~/bin/trilinos-12.10.1/" ..
