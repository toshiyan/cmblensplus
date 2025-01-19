#!/bin/bash

# define directory
dir_root=$(pwd)
dir_src_wrapped=${dir_root}/fortran_wrapped
dir_script=${dir_root}/scripts

# define functions
interface()
{
  echo 'creating interface python scripts with docstring' ${1}
  cd ${dir_src_wrapped}/src_${1}/
  python ${dir_script}/generate_docstring.py -libname lib${1} -modname ${2}
  mv *.py ${dir_root}/cmblensplus/${1}/
  cd ${dir_root}
}

# create python interface
if [ ${1} = "basic" -o ${1} = "all" ]; then
  interface basic "wigner_funcs.f90 bispec.f90 flat.f90 delens.f90 galaxy.f90 cosmofuncs.f90"
fi
if [ ${1} = "curvedsky" -o ${1} = "all" ]; then
  interface curvedsky "utils.f90 cninv.f90 norm_quad.f90 norm_imag.f90 rec_lens.f90 rec_ilens.f90 rec_rot.f90 rec_tau.f90 rec_src.f90 rec_iamp.f90 delens.f90 bispec.f90"
fi
if [ ${1} = "flatsky" -o ${1} = "all" ]; then
  interface flatsky "utils.f90 bispec.f90 rec_lens.f90 norm_lens.f90 rec_rot.f90 norm_rot.f90 rec_tau.f90 norm_tau.f90 rec_src.f90 norm_src.f90 norm_kxt.f90 norm_kxs.f90"
fi

