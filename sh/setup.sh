#!/bin/bash

fflags="-qopenmp -fPIC"
f2pycomp="intelem"

# f90 shared library
target=lib${pylib}

# set directory
flibloc="../../F90/"
lapack="${flibloc}/pub/LAPACK95/"
fftw="${flibloc}/pub/FFTW/"
cfitsio="${flibloc}/pub/cfitsio"
healpix="${flibloc}/pub/Healpix"
lenspix="${flibloc}/pub/lenspix"

pylibloc="../../py/"
wrapdir="../../wrap${pyv}/"

