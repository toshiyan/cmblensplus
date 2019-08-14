#!/bin/sh

pylib=basic

source ../../sh/pyversion.sh
source ../../sh/setup.sh

# compile link
modd="-I${flibloc}/mod"
libd="-L${flibloc}/lib"
link="-lutils"
# liomp5 and lpthread are required to use healpix
option="${modd} ${libd} ${link}"

# files to be compiled
scan="aps.f90 bispec.f90 flat.f90 delens.f90"

source ../../sh/compile.sh ${1}


