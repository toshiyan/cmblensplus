Installation Guide
==================

Requirements
------------

The package depends on Fortran Libraries, FFTW, HEALPix, cfitsio and LAPACK. 
Source code for those is included with this package. You also need F2py and a fortran compiler to compile Fortran sources.
Note that the F90 sources are written assuming the intel Fortran compiler. 

Install
-------

0) The pakage can be downloaded from https://github.com/toshiyan/cmblensplus or type 

    git clone git@github.com:toshiyan/cmblensplus.git


1) Edit the shell file to add, e.g. for .bashrc.ext in NERSC

    module load intel

    export LD_LIBRARY_PATH=${path-to-cmblensplus}/fortran_internal/src_public/Healpix/lib/:${path-to-cmblensplus}/fortran_internal/src_public/cfitsio/lib/:$LD_LIBRARY_PATH
    
    export PYTHONPATH=${path-to-cmblensplus}/utils/:${path-to-cmblensplus}/:$PYTHONPATH


2) python setup.py build

This automatically install the public codes (FFTW, cfitsio, Healpix, and Lapack), and then the local sources. 
Then the script produces python modules using f2py. 
The python modules can be found in cmblensplus/.


Tips
----

- libifport.so.5 not found

  Indicating that ifort is not loaded appropriately. 

- Segmentation fault after installing everything

  Some Healpix subroutine causes a problem of stack memory size. One quick solution is to set "ulimit -s unlimited"



