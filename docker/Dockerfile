# Get firedrake and install in same image as dolfin master
#
# Author:
# Chris Richardson
FROM quay.io/fenicsproject/stable:latest
MAINTAINER Chris Richardson <chris@bpi.cam.ac.uk>

USER root
WORKDIR /tmp
RUN apt-get update && \
    apt-get install -y automake git-core libopenmpi-dev libtool mercurial openmpi-bin python-pip python-tk && \
    pip2 install virtualenv

USER fenics
WORKDIR /home/fenics
RUN curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install && \
    unset PETSC_DIR && \
    unset PYTHONPATH && \
    python firedrake-install --disable-ssh
    
USER root
