#!/bin/sh

# setup directory

cwd=$(pwd)

fftw=${cwd}/FFTW
cfitsio=${cwd}/cfitsio
healpix=${cwd}/Healpix
lenspix=${cwd}/lenspix
lapack=${cwd}/LAPACK95


for args in "$@"
do 

  # FFTW
  if [ ${args} = "FFTW" -o ${args} = "all" ]; then
    echo '---- Install FFTW ----'
    cd ${fftw}
    wget -d https://www.fftw.org/fftw-3.3.8.tar.gz
    gunzip fftw-3.3.8.tar.gz
    tar -xvf fftw-3.3.8.tar
    cd fftw-3.3.8
    ./configure --prefix=${fftw}/ --enable-openmp --enable-shared --enable-static
    make
    make install
    cd ../
    rm -rf fftw-3.3.8 fftw-3.3.8.tar
    cd ${cwd}
  fi

  # cfitsio
  if [ ${args} = "cfitsio" -o ${args} = "all" ]; then
    echo '---- Install cfitsio ----'
    cd ${cfitsio}
    rm -rf bin lib include
    wget -d https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.5.0.tar.gz
    tar -xzf cfitsio-4.5.0.tar.gz
    cd cfitsio-4.5.0
    ./configure --prefix=${cfitsio}/
    make 
    make install
    cd ../
    rm -rf cfitsio-4.5.0 cfitsio-4.5.0.tar.gz
    cd ${cwd}
  fi

  # Healpix
  if [ ${args} = "healpix" -o ${args} = "all" ]; then
    echo '---- Install healpix ----'
    # Healpix (here, version 3.80)
    # In this version, -lsharp option is needed for compiling external source codes
    rm -rf Healpix
    mkdir Healpix
    wget -d https://sourceforge.net/projects/healpix/files/Healpix_3.80/Healpix_3.80_2021Jun22.tar.gz 
    tar xf Healpix_3.80_2021Jun22.tar.gz -C Healpix --strip-components 1
    cd Healpix
    # set path
    healpix_path=${healpix}/lib/
    cfitsio_path=${cfitsio}/lib/
    export LD_LIBRARY_PATH=${healpix_path}:${cfitsio_path}:$LD_LIBRARY_PATH
    echo "please edit bashrc.ext etc to add the following for future use:"
    echo "export LD_LIBRARY_PATH=${healpix_path}:${cfitsio_path}:\$LD_LIBRARY_PATH"
    # clean duplicate path
    export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | awk '!seen[$0]++' | tr '\n' ':' | sed 's/:$//')
    # --- install main codes --- #
    FC=ifort CC=gcc FITSDIR=${cfitsio}/lib ./configure --auto=f90
    make
    make test
    cd ${cwd}
    rm -rf Healpix_3.80_2021Jun22.tar.gz
  fi

  # lenspix
  if [ ${args} = "lenspix" -o ${args} = "all" ]; then
    echo '---- Install lenspix ----'
    cd ${lenspix}/src/
    make clean; make; make install
    cd ${cwd}
  fi

  # LAPACK (with ifort)
  if [ ${args} = "lapack" -o ${args} = "all" ]; then
    echo '---- Install LAPACK and LAPACK95 ----'
    cd ${lapack}
    # first install LAPACK
    mkdir tmp0 tmp1
    wget -d http://www.netlib.org/lapack/lapack.tgz ; tar zxvf lapack.tgz -C tmp0 --strip-components 1
    wget -d http://www.netlib.org/lapack95/lapack95.tgz ; tar zxvf lapack95.tgz -C tmp1 --strip-components 1
    rm -rf lapack.tgz lapack95.tgz
    cd ${lapack}/tmp0/
    cp INSTALL/make.inc.ifort make.inc
    sed -i "s/FFLAGS =/FFLAGS = -fPIC/g" make.inc
    sed -i "s/FFLAGS_NOOPT =/FFLAGS_NOOPT = -fPIC/g" make.inc
    make
    mv *.a ${lapack}/lib/
    # next install LAPACK95
    cd ${lapack}/tmp1/
    sed -i "s/ = f95 -free/ = ifort -fPIC/g" make.inc
    sed -i "s/ = f95 -fixed/ = ifort -fPIC/g" make.inc
    cd ${lapack}/tmp1/SRC/
    make single_double_complex_dcomplex
    mv ${lapack}/tmp1/lapack95.a ${lapack}/lib/liblapack95.a
    mv ${lapack}/tmp1/lapack95_modules/* ${lapack}/mod/lapack95_modules/
    # clean up
    cd ${lapack}
    rm -rf tmp0 tmp1
    cd ${cwd}
  fi

done

