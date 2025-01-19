#!/bin/bash

fflags="-qopenmp -fPIC"
#fflags="-qopenmp -fPIC -check -traceback"
f2pycomp="intelem"

# f90 shared library
target=lib${pylib}

#//// set directory ////#
# use installed public fortran modules in F90/pub/
cmblensplus="/global/homes/t/toshiyan/Work/Lib/cmblensplus/"

flibloc="${cmblensplus}F90/"
lapack="${flibloc}/pub/LAPACK95/"
fftw="${flibloc}/pub/FFTW/"
cfitsio="${flibloc}/pub/cfitsio"
healpix="${flibloc}/pub/Healpix"
lenspix="${flibloc}/pub/lenspix"
# use libraris already installed in the system
#fftw="/opt/cray/pe/fftw/3.3.10.6/x86_64/" # for NERSC cray-fftw

pylibloc="../../py/"
wrapdir="../../wrap${pyv}/"

