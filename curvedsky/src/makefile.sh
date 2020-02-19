#!/bin/sh

pylib=curvedsky

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${lenspix}/mod -I${healpix}/include -I${lapack}/mod -I${flibloc}/mod"
libd="-L${lenspix}/lib -L${healpix}/lib -L${lapack}/lib -L${cfitsio} -L${flibloc}/lib"
link="-lhp -llenspix -lmatrix -lutils -lhealpix -lcfitsio -liomp5 -lpthread -llapack95 -llapack -lrefblas"
# liomp5 and lpthread are required to use healpix

# files to be compiled
#scan="utils.f90 cninv.f90"
scan="utils.f90 cninv.f90 rec_lens.f90 norm_lens.f90 rec_rot.f90 norm_rot.f90 rec_tau.f90 norm_tau.f90 rec_src.f90 norm_src.f90 delens.f90 bispec.f90"

source ../../sh/compile.sh ${1}

