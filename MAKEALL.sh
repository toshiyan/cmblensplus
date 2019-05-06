#!/bin/sh

cwd=$(pwd)

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

if [ ${1} = "F90" -o ${1} = "all" ]; then
  rm -rf F90/lib/*.a
  rm -rf F90/mod/*.mod
  makeF90 src_utils
  makeF90 src_matrix
  makeF90 src_dft
fi

if [ ${1} = "f2py" -o ${1} = "all" ]; then
  #rm -rf wrap/*.so
  #rm -rf wrap_py2/*.so
  makef2py basic
  makef2py flatsky
  makef2py curvedsky
fi


