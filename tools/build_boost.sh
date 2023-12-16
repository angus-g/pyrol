#!/bin/sh

VERSION=1.83.0

wget "https://github.com/boostorg/boost/releases/download/boost-${VERSION}/boost-${VERSION}.tar.xz"
tar xf boost-${VERSION}.tar.xz
cd boost-${VERSION}
echo "using darwin : : $CXX ;" > user-config.jam
./bootstrap.sh --prefix=/usr/local
./b2 --prefix=/usr/local --user-config=user-config.jam