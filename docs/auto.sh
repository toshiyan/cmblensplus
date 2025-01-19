#!/bin/bash
rm -rf html
make html
mv build/html/ .

rename _s s html/*

find html/ -name "*.html" | xargs sed -i -e "s/_static/static/g" ${1}
find html/ -name "*.html" | xargs sed -i -e "s/_sources/sources/g" ${1}

git add *
git commit -m "testing docstring"
git push

