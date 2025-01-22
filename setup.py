from setuptools import Command, setup, Extension
#from setuptools.command.build_py import build_py as _build
#from setuptools.command.build_ext import build_ext
from numpy.distutils.core import setup, Extension # python v3.12 and later does not have this function...
from numpy.distutils.command.build_ext import build_ext as _build
#from setuptools.command.clean import clean
import numpy
import subprocess
import os
import glob
import shutil


os.environ['F90'] = 'ifort'
os.environ['F77'] = 'ifort'
os.environ['F90FLAGS'] = '-O3'  # Add optimization flags for Intel Fortran

# define module name
modname = ['basic','curvedsky','flatsky']


# Paths to directories that contain fortran sources
dir_fortran_int = 'fortran_internal/'
dir_src_utils   = dir_fortran_int + 'src_utils/'
dir_src_dft     = dir_fortran_int + 'src_dft/'
dir_src_hp      = dir_fortran_int + 'src_hp/'
dir_src_matrix  = dir_fortran_int + 'src_matrix/'

dir_src = {key: f"fortran_wrapped/src_{key}/" for key in modname}

dir_fortran_pub = dir_fortran_int + 'src_public/'
dir_fftw    = dir_fortran_pub + 'FFTW/'
dir_cfitsio = dir_fortran_pub + 'cfitsio/'
dir_lapack  = dir_fortran_pub + 'LAPACK95/'
dir_healpix = dir_fortran_pub + 'Healpix/'


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

# Define Extensions
extensions = [
    Extension(
        name='libbasic',
        sources=sources['basic'],
        libraries=["utils"],
        library_dirs=[dir_fortran_int+'/lib'],
        include_dirs=[dir_fortran_int+'/mod', numpy.get_include()],
        extra_compile_args=["-O3"],
    ),
    Extension(
        name='libcurvedsky',
        sources = sources['curvedsky'], # List of Fortran files in curvedsky/src
        libraries = ['hp','matrix','utils','healpix','sharp','cfitsio','iomp5','pthread','lapack95','lapack','refblas'],  
        library_dirs=[dir_fortran_int+'/lib',dir_healpix+'/lib',dir_cfitsio+'/lib',dir_lapack+'/lib'],
        include_dirs=[dir_fortran_int+'/mod',dir_healpix+'/include',dir_lapack+'/mod', numpy.get_include()],
        #extra_compile_args=["-g", "-O0"],
        extra_compile_args=["-O3"],
        extra_link_args=[""],  # Optional: Additional linker flags
    ),
    Extension(
        name='libflatsky',
        sources=sources['flatsky'],
        libraries = ['dft','utils','fftw3'],  
        library_dirs=[dir_fortran_int+'/lib',dir_fftw+'/lib'],
        include_dirs=[dir_fortran_int+'/mod',dir_fftw+'/include',numpy.get_include()],
        extra_compile_args=["-O3"],
    ),
]


#class CustomBuildExt(build_ext):
class CustomBuildExt(_build):
    
    def run(self):
        # Custom function to run before compiling
        print('run make')
        # compile all local f90 files in fortran_internal to create a local static library file
        #self.install_public_f90_library()
        self.install_local_f90_library()
        # Generate .f90 files from .src in ./fortran_wrapped/ with f2py directive
        self.add_f2py_directive_in_fortran_wrapped()
        # generate interface python files in ./cmblensplus/ for each submodule
        self.generate_interface_to_cmblensplus()
        # Call the original build_ext.run() to proceed with the actual build process
        print('compile sources inside fortran_wapped/')
        super().run()

    def build_extension(self, ext):
        ext_path = self.get_ext_fullpath(ext.name) # Ensure the .so file is saved in the package directory
        ext_dir = os.path.dirname(ext_path)
        ext_dir = os.path.join(os.getcwd(), "cmblensplus")
        os.makedirs(ext_dir, exist_ok=True)
        self.build_lib = ext_dir
        ext._full_output_dir = ext_dir
        super().build_extension(ext)

        return CustomBuildExt

    def install_public_f90_library(self):
        subprocess.check_call(["cd",dir_fortran_pub])
        subprocess.check_call(["./install.sh", "FFTW, cfitsio, healpix, lapack"])
        subprocess.check_call(["cd","../../"])

    def install_local_f90_library(self):
        for directory in [dir_src_utils, dir_src_dft, dir_src_matrix, dir_src_hp]:
            subprocess.check_call(["make", "-C", directory])
            subprocess.check_call(["make", "install", "-C", directory])

    def add_f2py_directive_in_fortran_wrapped(self):
        for mod in modname:
            for srcname in sources[mod]:
                subprocess.check_call(["python", "scripts/f2pysig.py", "-srcname", srcname.replace('.f90','.src')])

    def generate_interface_to_cmblensplus(self):
        subprocess.check_call(["sh", "scripts/create_interface.sh", "all"])
    

# Custom clean command
class CustomClean(Command):

    user_options = []  # Required attribute: empty list for no options

    def initialize_options(self):
        """Set default values for options (required)."""
        pass

    def finalize_options(self):
        """Post-process options (required)."""
        pass

    def run(self):
        paths_to_clean = ["build", "dist", "*.egg-info", "fortran_internal/**/*.o"]
        for path in paths_to_clean:
            for file_or_dir in glob.glob(path, recursive=True):
                try:
                    if os.path.isfile(file_or_dir):
                        os.remove(file_or_dir)
                    elif os.path.isdir(file_or_dir):
                        shutil.rmtree(file_or_dir)
                    print(f"Removed: {file_or_dir}")
                except Exception as e:
                    print(f"Error removing {file_or_dir}: {e}")

# Setup configuration
setup(
    name = 'cmblensplus',
    ext_modules=extensions,
    include_dirs = [numpy.get_include()],
    cmdclass = {
        'build_ext': CustomBuildExt,
        'clean': CustomClean,
    },
)

