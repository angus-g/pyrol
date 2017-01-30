INCLUDES=-I/usr/include/python2.7 -I$(TRILINOS_DIR)/include -I/usr/include/mpi
OPTS=-std=c++11 -shared -fPIC
LIBS=-lmpi -lmpi_cxx -lrol -lteuchoscore -lteuchoscomm -lpython2.7 -L$(TRILINOS_DIR)/lib
all : ROL.so
ROL.so : ROL_wrap.cxx
	g++ $(OPTS) -o ROL.so ROL_wrap.cxx $(INCLUDES) $(LIBS)
