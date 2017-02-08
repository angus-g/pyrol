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
    curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install && \
    sudo -H python3 firedrake-install
