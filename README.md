# Analysis Tools for 2D Maps and Higher-Order Statistics

`cmblensplus` is a Python package for reconstructing the lensing potential, cosmic birefringence, and patchy reionization from cosmic microwave background anisotropies (CMB), in both full-sky and flat-sky geometries. The package also includes tools for delensing, lensing bispectrum calculations, and related analyses.

Installation on NERSC systems is straightforward.

# Installation

The default installation builds the `basic` module and installs the pure-Python `curvedsky` and `utils` modules. This default installation does not require external Fortran libraries.

```bash
python -m pip install "cmblensplus @ git+https://github.com/toshiyan/cmblensplus.git@main"
```

For an editable install from a local checkout, run the following command from the top-level source directory:

```bash
python -m pip install --no-build-isolation --no-cache-dir -e ".[dev]"
```

## Installing the optional `flatsky` module

The `flatsky` module depends on FFTW. To install this module, first load the FFTW module in your environment. At NERSC, you can use:

```bash
module load cray-fftw
```

On NERSC systems, it is also recommended to use the Cray compiler wrappers when building the Fortran extensions:

```bash
export FC=ftn
export F90=ftn
export F77=ftn
export CC=cc
```

Then install the package with the optional `flatsky` module enabled:

```bash
python -m pip install -Csetup-args=-Dflatsky=true \
  "cmblensplus @ git+https://github.com/toshiyan/cmblensplus.git@main"
```

For an editable install with the optional `flatsky` module enabled, use:

```bash
python -m pip install --no-build-isolation --no-cache-dir -e ".[dev]" \
  -Csetup-args=-Dflatsky=true
```

Without `-Dflatsky=true`, the package is installed without the FFTW-dependent `flatsky` extension.

The package modules are located under `cmblensplus/`.

## Dependencies

The main dependencies are:

```text
numpy
scipy
matplotlib
astropy
healpy
ducc0
meson
ninja
meson-python
pkg-config
compilers
ipykernel
pip
tqdm
```

Some example Jupyter notebooks use `pytempura` to compute the quadratic estimator normalization:

https://github.com/simonsobs/tempura

## Build Process

The installation process consists of the following steps:

1. Compile local Fortran codes to create internal libraries under `fortran_internal/src_*/`.
2. Compile Fortran wrappers under `fortran_wrapped/src_*/` using `f2py`.
3. Build the Python extension modules under `cmblensplus/*/`.

# Documentation and References

Documentation is available at:

https://toshiyan.github.io/clpdoc/html/

The reference papers for each module are listed below.

## Curved-Sky Modules

* **Lensing Reconstruction and Delensing**
  Developed by Namikawa & Nagata, JCAP 09 (2014) 009
  https://arxiv.org/abs/1405.6568

* **Cosmic Birefringence**
  Developed by Namikawa et al., PRD 101 (2020) 083527
  https://arxiv.org/abs/2001.10465

* **Patchy Reionization**
  Developed by Namikawa, PRD 97 (2018) 063505
  https://arxiv.org/abs/1711.00058

* **Lensing Bispectrum**
  Developed by Namikawa et al., PRD 99 (2019) 063511
  https://arxiv.org/abs/1812.10635

## Flat-Sky Modules

* **Lensing Reconstruction and Delensing**
  Developed by Namikawa, PRD 95 (2017) 103514
  https://arxiv.org/abs/1703.00169

* **Cosmic Birefringence**
  Developed by Namikawa, PRD 95 (2017) 043523
  https://arxiv.org/abs/1612.07855

# Package Structure

This package contains four main Python modules based on Fortran 90 source codes and pure-python codes:

* `cmblensplus.basic`
  Basic routines, such as analytic calculations of delensed B-mode spectra and lensing bispectra.

* `cmblensplus.curvedsky`
  Tools for full-sky analyses, including lensing reconstruction, cosmic birefringence, patchy reionization optical depth, bias hardening, bispectrum calculations, delensing, and analytic reconstruction normalization. The normalization routines are no longer actively developed within this software and are now developed in `pytempura`.

* `cmblensplus.flatsky`
  Tools for flat-sky analyses corresponding to the curved-sky modules. This module requires FFTW.

Additional Python scripts are provided in:

* `cmblensplus.utils`
  Utility scripts to compute noise bias, RDN0, diagonal RDN0, and related quantities.


# Examples

Example scripts are available in the `example/` directory.

# Acknowledgements

This package uses the external FFTW library for the `flatsky` module.

# Contact

Toshiya Namikawa
toshiya.namikawa at ipmu.jp
