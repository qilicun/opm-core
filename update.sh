#!/bin/sh
echo "starting fetch OPM"
git fetch opm
git checkout master
git merge opm/master
git push
echo "Starting build OPM"
cd build/
nice make -j 4
make tests
make exmaples
make doc
make install
echo "All things going fine!"
