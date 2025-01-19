from distutils.errors import DistutilsError
from numpy.distutils.core import setup, Extension # python v3.12 and later does not have this function...
from numpy.distutils.command.build_ext import build_ext as numpy_build_ext
from setuptools import Command
import numpy
import subprocess
import os
import glob


# define module name
modname = ['basic','curvedsky','flatsky']


# Paths to directories that contain fortran sources
dir_fortran_int = 'fortran_internal/'
dir_src_utils   = dir_fortran_int + 'src_utils/'
dir_src_dft     = dir_fortran_int + 'src_dft/'
dir_src_hp      = dir_fortran_int + 'src_hp/'
dir_src_matrix  = dir_fortran_int + 'src_matrix/'
lib_fortran_int = dir_fortran_int + 'lib/'
mod_fortran_int = dir_fortran_int + 'mod/'

dir_fortran_wrp = 'fortran_wrapped/'
dir_src = {key: f"{dir_fortran_wrp}src_{key}/" for key in modname}

dir_fortran_pub = dir_fortran_int + 'src_public/'
dir_fftw    = dir_fortran_pub + 'FFTW/'
dir_cfitsio = dir_fortran_pub + 'cfitsio/'
dir_lapack  = dir_fortran_pub + 'LAPACK95/'
dir_healpix = dir_fortran_pub + 'Healpix/'
dir_lenspix = dir_fortran_pub + 'lenspix/'


# define source files of each module
fname = {
    'basic':     ['wigner_funcs','bispec','flat','delens','galaxy','cosmofuncs'],
    'curvedsky': ['utils','cninv','norm_quad','norm_imag','rec_lens','rec_ilens','rec_rot','rec_tau','rec_src','rec_iamp','delens','bispec'], 
    'flatsky':   ['utils','bispec','rec_lens','norm_lens','rec_rot','norm_rot','rec_tau','norm_tau','rec_src','norm_src','norm_kxt','norm_kxs']
}


# source files
sources = {
    mod: [ dir_src[mod] + f'{srcname}.f90' for srcname in fname[mod] ] for mod in modname
}

print(sources['basic'])

# Function to run Makefile for compiling internal fortran codes and for generating static libraries
def run_make():
    """Run the Makefile to compile utils.a and hp.a."""
    dirs = [dir_src_utils, dir_src_dft, dir_src_matrix, dir_src_hp]
    for directory in dirs:
        subprocess.check_call(["make", "-C", directory])
        subprocess.check_call(["make", "install", "-C", directory])


def add_f2py_directive():
    for mod in modname:
        for srcname in sources[mod]:
            subprocess.check_call(["python", "scripts/f2pysig.py", "-srcname", srcname.replace('.f90','.src')])

def generate_interface():
    subprocess.check_call(["sh", "scripts/create_interface.sh", "all"])


# Run Makefile to build the static libraries before continuing with setup.py
run_make()

# Generate .f90 files with f2py directive
add_f2py_directive()

# generate interface python files for each submodule
generate_interface()

# Define the Fortran extensions
libbasic_extension = Extension(
    name = 'libbasic',  # The Python module name
    sources = sources['basic'], # List of Fortran files in curvedsky/src
    libraries=["utils"],  # Link to the static libraries (utils.a and hp.a)
    library_dirs=[lib_fortran_int],
    include_dirs=[mod_fortran_int, numpy.get_include()],
    extra_compile_args=["-O3"],  # Optional: Optimization flags for the Fortran compiler
    extra_link_args=[],  # Optional: Additional linker flags
)

libcurvedsky_extension = Extension(
    name = 'libcurvedsky',  # The Python module name
    sources = sources['curvedsky'], # List of Fortran files in curvedsky/src
    libraries = ['hp','matrix','utils','lenspix','healpix','sharp','cfitsio','iomp5','pthread','lapack95','lapack','refblas'],  
    library_dirs=[lib_fortran_int,dir_healpix+'/lib',dir_cfitsio+'/lib',dir_lapack+'/lib',dir_lenspix+'/lib'],
    include_dirs=[mod_fortran_int,dir_healpix+'/include',dir_lapack+'/mod',dir_lenspix+'/mod', numpy.get_include()],
    #extra_compile_args=["-g", "-O0"],
    extra_compile_args=["-O3"],
    extra_link_args=[""],  # Optional: Additional linker flags
)

libflatsky_extension = Extension(
    name = 'libflatsky',  # The Python module name
    sources = sources['flatsky'], # List of Fortran files in curvedsky/src
    libraries = ['dft','utils','fftw3'],  
    library_dirs=[lib_fortran_int,dir_fftw+'/lib'],
    include_dirs=[mod_fortran_int,dir_fftw+'/include',numpy.get_include()],
    extra_compile_args=["-O3"],
    extra_link_args=[""],  # Optional: Additional linker flags
)


class CustomBuildExt(numpy_build_ext):
    
    def build_extension(self, ext):
        
        # Ensure the .so file is saved in the package directory
        ext_path = self.get_ext_fullpath(ext.name)
        ext_dir = os.path.dirname(ext_path)
        ext_dir = os.path.join(os.getcwd(), "cmblensplus")
        os.makedirs(ext_dir, exist_ok=True)
        self.build_lib = ext_dir
        ext._full_output_dir = ext_dir
        super().build_extension(ext)

        return CustomBuildExt


# Setup configuration
setup(
    name = 'cmblensplus',
    ext_modules = [libbasic_extension,libcurvedsky_extension,libflatsky_extension], # Register the extension
    include_dirs = [numpy.get_include()],
    cmdclass = {
        'build_ext': CustomBuildExt,  # Use CustomBuildExt to build extensions and move .so files
    },
)

