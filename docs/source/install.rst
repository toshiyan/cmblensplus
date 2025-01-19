Installation Guide
==================

Requirements
------------

The package depends on Fortran Libraries, FFTW, HEALPix, cfitsio and LensPix. 
Source code for those is included with this package. You also need F2py and a fortran compiler to compile Fortran sources.
Note that the F90 sources are written assuming the intel Fortran compiler. 

Install
-------

0) The pakage can be downloaded from https://github.com/toshiyan/cmblensplus or type 

    git clone git@github.com:toshiyan/cmblensplus.git


1) Edit the shell file to add, e.g. for .bashrc.ext in NERSC

    module load intel

    export LD_LIBRARY_PATH=${path-to-cmblensplus}/F90/pub/Healpix/lib/:${path-to-cmblensplus}/F90/pub/cfitsio/:$LD_LIBRARY_PATH
    
    export PYTHONPATH=${path-to-cmblensplus}/utils/:${path-to-cmblensplus}/wrap/:$PYTHONPATH


2) ./install.sh all 

The install.sh script automatically install the public codes (FFTW, cfitsio, Healpix, Lenspix, and Lapack) at F90/pub/, and then the local sources at F90/. 
Then the script produces python modules using f2py. 
The python modules can be found in wrap/ for python 3.
(wrap_py2 for python 2 is no longer supported)



Tips
----

- libifport.so.5 not found

  Indicating that ifort is not loaded appropriately. 

- Segmentation fault after installing everything

  Some Healpix subroutine has a problem. One quick solution is to set "ulimit -s unlimited"



