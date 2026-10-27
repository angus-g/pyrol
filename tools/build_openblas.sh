#!/bin/sh

VERSION=0.3.34

wget "https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${VERSION}.tar.gz"
tar xf v${VERSION}.tar.gz
cd OpenBLAS-${VERSION}

make CC=cc FC=gfortran-14 libs netlib shared
sudo make PREFIX=/usr/local install
