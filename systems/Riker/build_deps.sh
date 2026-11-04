#!/bin/bash
# c-lime and mpfr (no modules on Riker); hdf5 with C++ API only if the module lacks it.
# Run on a LOGIN node (downloads). Usage: source sourceme.sh; ./build_deps.sh
set -eo pipefail
mkdir -p $RIKER_SW/src && cd $RIKER_SW/src
if [ ! -f $DEPS/lib/liblime.a ]; then
  [ -d c-lime ] || git clone https://github.com/usqcd-software/c-lime.git
  (cd c-lime && ./autogen.sh && ./configure --prefix=$DEPS CC=gcc && make -j16 && make install)
fi
if [ ! -f $DEPS/include/mpfr.h ]; then
  [ -f mpfr-4.2.1.tar.xz ] || curl -LO https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz
  tar xf mpfr-4.2.1.tar.xz
  (cd mpfr-4.2.1 && ./configure --prefix=$DEPS --with-gmp=$GMP && make -j16 && make install)
fi
if [ ! -f $HDF5/include/H5Cpp.h ]; then
  echo "hdf5 module has no C++ API: building hdf5 1.14.6 --enable-cxx into $DEPS"
  [ -f hdf5-1.14.6.tar.gz ] || curl -LO https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
  tar xf hdf5-1.14.6.tar.gz
  (cd hdf5-1.14.6 && ./configure --prefix=$DEPS --enable-cxx CC=gcc CXX=g++ && make -j16 && make install)
fi
ls $DEPS/lib/liblime.a $DEPS/include/mpfr.h $HDF5/include/H5Cpp.h
