#!/bin/sh

# setup directory

cwd=$(pwd)
pub=${cwd}/F90/pub

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
    cd ${cfitsio}/cfitsio/
    make clean
    ./configure
    make 
    mv libcfitsio.a ../
    make clean
    cd ${cwd}
  fi

  # Healpix
  if [ ${args} = "healpix" -o ${args} = "all" ]; then
    echo '---- Install healpix ----'
    #cd ${healpix}
    #make clean
    #printf '%s \n %s \n %s \n %s \n %s \n %s \n %s \n %s \n %s \n %s \n %s' 3 ifort "" "" "" "" "" "" "" "" "../cfitsio"  | ./configure
    #make; make test
    #cd ${cwd}
    # Healpix (here, version 3.80)
    # from this version, -lsharp option is needed for compiling external source codes
    wget -d https://sourceforge.net/projects/healpix/files/Healpix_3.80/Healpix_3.80_2021Jun22.tar.gz 
    rm -rf Healpix
    mkdir Healpix
    tar xf Healpix_3.80_2021Jun22.tar.gz -C Healpix --strip-components 1
    cd Healpix
    ## --- install libsharp --- #
    ##edit Healpix/src/common_libraries/libsharp/configure as
    # --- install main codes --- #
    FC=ifort CC=icc FITSIO=../cfitsio ./configure --auto=f90
    make
    make test
    rm -rf Healpix_3.80_2021Jun22.tar.gz
    # --- edit bashrc.ext (the following is an example)
    # export LD_LIBRARY_PATH=$HOME/Work/Lib/cmblensplus/F90/pub/Healpix/lib/:$LD_LIBRARY_PATH

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

