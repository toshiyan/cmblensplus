# Analysis tools for 2D map to compute some higher-order statistics

This package contains a Fortran wrapper for Python, to reconstruct lensing potential, cosmic bi-refrimgence, and patchy reionization from cosmic microwave background anisotropies (CMB) in full and flat sky. This package also includes modules for delensing, some bi-spectrum calculation, and so on. Installing this library at NERSC is straightforward. 

# Documents and Reference

Please go to the following link for details of how to use each module:
https://toshiyan.github.io/clpdoc/html/. 
The reference paper for each module is as follows:

Curved sky modules:

  - Temperature / Polarization Lensing Reconstruction, Delensing
   developed in Namikawa & Nagata JCAP 09 (2014) 009, https://arxiv.org/abs/1405.6568
  - https://arxiv.org/abs/2001.10465 : Cosmic Birefringence
  - https://arxiv.org/abs/1711.00058 : Patchy Reionization
  - https://arxiv.org/abs/1812.10635 : Lensing bi-spectrum

Flat sky modules:

  - https://arxiv.org/abs/1209.0091 : Temperature Lensing Reconstruction, Bias-Hardening
  - https://arxiv.org/abs/1310.2372 : Temperature / Polarization Lensing Reconstruction
  - https://arxiv.org/abs/1612.07855 : Cosmic Bi-refringence
  - https://arxiv.org/abs/1703.00169 : CMB Delensing
  - https://arxiv.org/abs/1706.05133 : Lensing bi-spectrum


# This Package

This package contains three main python modules based on Fourtran 90 sources: 
  
  - basic     --- basic routines such as analytic calculation of delensed B-mode spectrum and lensing bispectrum.

  - curvedsky --- analysis package to measure lensing, birefringence, patchy tau, bias-hardening, bispectrum, delensing and analytic reconstruction normalization.
  
  - flatsky   --- the same as curvedsky code but for flatsky analysis.

The additional simple python sctipts are stored at

  - example   --- several example scripts for demonstration
  
  - utils     --- python scripts to compute noise biase, RDN0, diagonal RDN0, etc. 


# Installation

  [1] Compile Fortran public software

  Go to "F90/pub/" directory, and install each pubclic package (FFTW, Healpix, LAPACK, cfitsio and Lenspix). 

  [2] Create Fortran wrapper

  Go to the root directory and type one of the following command:

    - ./MAKEALL.sh all       (for generating all modules)
    - ./MAKEALL.sh basic     (for basic module)
    - ./MAKEALL.sh curvedsky (for curvedsky module)
    - ./MAKEALL.sh flatsky   (for flatsky module)
  
  You will find modules under "wrap/" for python 3 and "wrap_py2" for python 2, depending on your f2py version.

# Examples

You can find the example codes at "example" directory. 


# Acknowledgement

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO, LensPix. 

# Contact

  namikawa at slac.stanford.edu

