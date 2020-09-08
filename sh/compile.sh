#!/bin/bash

option="${modd} ${libd} ${link}"
obj="${target}.pyf ${scan}"

if [ ${1} != "all" -a ${1} != "scan" -a ${1} != "compile" -a ${1} != "docs" -a ${1} != "clean" -a ${1} != "dev" ]; then
	echo "the argument must be either all, scan, compile, docs, clean or dev"
fi

# create signature file
if [ ${1} = "scan" -o ${1} = "all" -o ${1} = "dev" ]; then
  rm -rf ${target}.pyf
  python ${pylibloc}/f2pysig.py -libname ${target} -modname ${scan}
fi 

# compile
if [ ${1} = "compile" -o ${1} = "all" -o ${1} = "dev" ]; then
  f2py --f90flags="${fflags}" --fcompiler=${f2pycomp} --noopt ${option} -c ${obj} -m ${target} 
  # save lib
  mv ${target}*.so ${wrapdir}/${target}.so
fi

# create python interface
if [ ${1} = "docs" -o ${1} = "all" ]; then
  python ${pylibloc}/docstr.py -libname ${target} -modname ${scan} ${2}
  mv *.py ${wrapdir}/${pylib}/
fi

if [ ${1} = "clean" -o ${1} = "all" ]; then
  rm -rf ${target}.pyf
fi


