#!/bin/sh
git fetch opm
git checkout master
git merge opm/master
git push

cd build/
nice make -j 4
make tests
make exmaples
make doc
make install
