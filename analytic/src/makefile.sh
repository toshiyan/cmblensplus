#!/bin/sh

pylib=analytic

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${lapack}/mod -I${flibloc}/mod"
libd="-L${lapack}/lib -L${flibloc}/lib"
link="-lmatrix -lutils -llapack95 -llapack -lrefblas"

# files to be compiled
scan="rec_lens.f90 rec_rot.f90 rec_tau.f90 rec_src.f90 delens.f90 flat.f90"

source ../../sh/compile.sh ${1}

