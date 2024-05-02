#!/bin/sh

VERSION=1.85.0
SUFFIX=-b2-nodocs

[ ! -f "boost-${VERSION}.tar.xz" ] && wget "https://github.com/boostorg/boost/releases/download/boost-${VERSION}/boost-${VERSION}${SUFFIX}.tar.xz"
[ -d "boost-${VERSION}" ] && rm -rf "boost-${VERSION}"
tar xf boost-${VERSION}${SUFFIX}.tar.xz
cd boost-${VERSION}
echo "using darwin : : $CXX ;" > user-config.jam
python_include=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["include"])')
echo "using python : : : $python_include ;" >> user-config.jam
./bootstrap.sh --prefix=/usr/local
sudo ./b2 --prefix=/usr/local --user-config=user-config.jam install
