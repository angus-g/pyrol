# Install DOLFIN and ROL/Trilinos together
#
# Author:
# Chris Richardson
FROM quay.io/fenicsproject/dev-env:latest
MAINTAINER Chris Richardson <chris@bpi.cam.ac.uk>

USER root
RUN pip3 install git+https://bitbucket.org/fenics-project/fiat.git@master && \
    pip3 install git+https://bitbucket.org/fenics-project/ufl.git@master && \
    pip3 install git+https://bitbucket.org/fenics-project/dijitso.git@master && \
    pip3 install git+https://bitbucket.org/fenics-project/ffc.git@master && \
    git clone https://bitbucket.org/fenics-project/dolfin.git && \
    cd dolfin && mkdir -p build && cd build && cmake .. && make -j2 && make install

RUN apt-get -qq update && \
    apt-get -y --with-new-pkgs \
        -o Dpkg::Options::="--force-confold" upgrade && \
    apt-get -y install patchelf
RUN cd dolfin/python && pip3 install . && pip3 install roltrilinos
RUN unset PETSC_DIR; curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install; python3 firedrake-install

