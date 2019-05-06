# CMB data analysis tool for some higher-order statistics

This package contains a Fortran wrapper for Python, to reconstruct lensing potential, cosmic bi-refrimgence, and patchy reionization from cosmic microwave background anisotropies (CMB) in full and flat sky. This package also includes modules for delensing, some bi-spectrum calculation, and so on. Installing this library at NERSC is straightforward. 

# Documents

Please go to the following link for details:
https://toshiyan.github.io/clpdoc/html/

# Installation

  [1] Compile Fortran public software

  Go to "F90/pub/" directory, and install each pubclic package (FFTW, Healpix, LAPACK, cfitsio and Lenspix). 

  [2] Create Fortran wrapper

  Go to the root directory and type ./"MAKEALL.sh all". You will find modules under "wrap/" for python 3 and "wrap_py2" for python 2, depending on your f2py version. 

# Examples

You can find the example codes at "example" directory. 


# References

In flat sky, the source code has been used for

  - https://arxiv.org/abs/1209.0091 : the temperature lensing reconstruction, 
  - https://arxiv.org/abs/1310.2372 : the polarization lensing reconstruction, 
  - https://arxiv.org/abs/1612.07855 : the cosmic bi-refringence reconstruction, 
  - https://arxiv.org/abs/1703.00169 : the CMB delensing
  - https://arxiv.org/abs/1706.05133 : the bi-spectrum

The curvedsky analysis tools are also tested as follows for the lensing reconstruction, patchy reionization reconstruction, bi-spectrum and delensing

  - https://arxiv.org/abs/1405.6568 : the temperature/polarization lensing reconstruction, and delensing
  - https://arxiv.org/abs/1711.00058 : the patchy reionization reconstruction
  - https://arxiv.org/abs/1812.10635 : the bi-spectrum


# Acknowledgement

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO, LensPix. 

# Contact

  namikawa at slac.stanford.edu

