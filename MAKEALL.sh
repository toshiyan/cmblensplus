#!/bin/sh

cwd=$(pwd)

# define functions
makeF90()
{
  echo '----' ${1} '----'
  cd F90/${1}
  make -f Makefile clean
  make -f Makefile
  make -f Makefile install
  make -f Makefile clean
  cd ${cwd}
}

makef2py()
{
  echo '----' ${1} '----'
  cd ${1}/src/
  ./makefile.sh all
  cd ${cwd}
}

if [ ${1} = "clean" ]; then
  rm -rf wrap/*.so
  rm -rf wrap_py2/*.so
  rm -rf F90/lib/*.a
  rm -rf F90/mod/*.mod
fi

# compile f90 sources
if [ ${1} = "F90" -o ${1} = "all" -o ${1} = "basic" -o ${1} = "flatsky" -o ${1} = "curvedsky" ]; then
  rm -rf F90/lib/*.a
  rm -rf F90/mod/*.mod
  makeF90 src_utils
  makeF90 src_matrix
fi
if [ ${1} = "F90" -o ${1} = "all" -o ${1} = "flatsky" ]; then
  makeF90 src_dft
fi
if [ ${1} = "F90" -o ${1} = "all" -o ${1} = "curvedsky" ]; then
  makeF90 src_hp
fi

# create python modules
if [ ${1} = "f2py" -o ${1} = "all" -o ${1} = "basic" ]; then
  makef2py basic
fi
if [ ${1} = "f2py" -o ${1} = "all" -o ${1} = "flatsky" ]; then
  makef2py flatsky
fi
if [ ${1} = "f2py" -o ${1} = "all" -o ${1} = "curvedsky" ]; then
  makef2py curvedsky
fi

