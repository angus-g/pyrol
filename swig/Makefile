INCLUDES=-I/usr/include/python2.7 -I$(TRILINOS_DIR)/include -I/usr/include/mpi
OPTS=-std=c++11 -shared -fPIC
LIBS=-lmpi -lmpi_cxx -lrol -lteuchoscore -lteuchoscomm -lteuchosnumerics -lpython2.7 -L$(TRILINOS_DIR)/lib
all : _ROL.so
ROL_wrap.cxx : ROL.i Makefile ROLVector.h ROLObjective.h
	swig $(INCLUDES) -c++ -python $<

_ROL.so : ROL_wrap.cxx
	g++ $(OPTS) -o $@ ROL_wrap.cxx $(INCLUDES) $(LIBS)
