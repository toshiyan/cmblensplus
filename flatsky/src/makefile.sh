#!/bin/sh

pylib=flatsky

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${fftw}/api -I${flibloc}/mod"
libd="-L${fftw} -L${flibloc}/lib"
link="-ldft -lutils -lfftw3"

# files to be compiled
#scan="rec_lens.f90 norm_lens.f90 rec_rot.f90 rec_tau.f90 delens.f90 remap.f90 bispec.f90"
scan="utils.f90 ffttools.f90 rec_lens.f90 norm_lens.f90"

source ../../sh/compile.sh ${1}

