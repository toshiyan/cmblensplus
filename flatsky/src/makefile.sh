#!/bin/sh

pylib=flatsky

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${fftw}/api -I${flibloc}/mod"
libd="-L${fftw} -L${flibloc}/lib"
link="-ldft -lutils -lfftw3"

# files to be compiled
scan="utils.f90 ffttools.f90 rec_lens.f90 norm_lens.f90 rec_rot.f90 norm_rot.f90 rec_tau.f90 norm_tau.f90 rec_src.f90"

source ../../sh/compile.sh ${1}

