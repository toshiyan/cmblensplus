#!/bin/sh

version=$(python -V 2>&1 | grep -Po '(?<=Python )(.+)')
if [[ "${version::1}" -lt "3" && "${version::1}" -gt "1" ]]

then
  pyv=_py2
  echo $pyv
else
  pyv=""
fi


