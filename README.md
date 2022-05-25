# Analysis tools for 2D map to compute some higher-order statistics

This package contains a wrapper for Python, to reconstruct lensing potential, cosmic bi-refrimgence, and patchy reionization from cosmic microwave background anisotropies (CMB) in full and flat sky. This package also includes modules for delensing, some bi-spectrum calculation, and so on. Installing this library at NERSC is straightforward. 

# Installation

  The easiest way to install the entire package is to run the shellscript: 
     - ./install.sh all
  You will find modules inside "wrap/". 

  Note that the install.sh file compiles the following files:  

  [1] Fortran public software (located inside F90/pub/)

  FFTW, cfitsio, Healpix, LAPACK95, and Lenspix. 

  [2] Fortran codes to create a wrapper


# Documents and Reference

Please go to the following link for details of how to use each module:
https://toshiyan.github.io/clpdoc/html/. 
The reference paper for each module is as follows:

**Curved sky modules**:

  - **Lensing Reconstruction, Delensing** \
   Developed by Namikawa & Nagata JCAP 09 (2014) 009, https://arxiv.org/abs/1405.6568
  - **Cosmic Birefringence** \
   Developed by Namikawa et al. PRD 101 (2020) 083527, https://arxiv.org/abs/2001.10465
  - **Patchy Reionization** \
   Developed by Namikawa PRD, 97 (2018) 063505, https://arxiv.org/abs/1711.00058
  - **Lensing Bi-Spectrum** \
   Developed by Namikawa et al. PRD 99 (2019) 063511, https://arxiv.org/abs/1812.10635

**Flat sky modules**:

  - **Lensing Reconstruction, Delensing** \
   Developed by Namikawa PRD 95 (2017) 103514, https://arxiv.org/abs/1703.00169
  - **Cosmic Birefringence** \
   Developed by Namikawa PRD 95 (2017) 043523, https://arxiv.org/abs/1612.07855


# This Package

This package contains three main python modules based on Fourtran 90 sources: 
  
  - basic     --- basic routines such as analytic calculation of delensed B-mode spectrum and lensing bispectrum.

  - curvedsky --- analysis package to measure lensing, birefringence, patchy tau, bias-hardening, bispectrum, delensing and analytic reconstruction normalization.
  
  - flatsky   --- the same as curvedsky code but for flatsky analysis.

The additional simple python sctipts are stored inside

  - example   --- several example scripts for demonstration
  
  - utils     --- python scripts to compute noise biase, RDN0, diagonal RDN0, etc. 


# Examples

You can find example codes inside "example" directory. 


# Acknowledgement

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO, LensPix. 

# Contact

  namikawa at slac.stanford.edu

