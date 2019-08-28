#!/bin/sh

pylib=curvedsky

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${lenspix}/mod -I${healpix}/include -I${flibloc}/mod"
libd="-L${lenspix}/lib -L${healpix}/lib -L${cfitsio} -L${flibloc}/lib"
link="-lhp -llenspix -lutils -lhealpix -lcfitsio -liomp5 -lpthread"
# liomp5 and lpthread are required to use healpix

# files to be compiled
scan="utils.f90 cninv.f90 rec_lens.f90 norm_lens.f90 rec_rot.f90 norm_rot.f90 rec_tau.f90 norm_tau.f90 rec_src.f90 norm_src.f90 delens.f90 bispec.f90"
#scan="utils.f90 rec_lens.f90 norm_lens.f90 rec_rot.f90 norm_rot.f90 rec_tau.f90 norm_tau.f90 rec_src.f90 norm_src.f90 delens.f90 bispec.f90"

source ../../sh/compile.sh ${1}

