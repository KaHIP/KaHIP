#!/bin/bash

rm -rf deploy
NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi
scons program=library variant=optimized -j $NCORES

mkdir deploy
cp ./optimized/interface/lib* deploy/
cp ./interface/kaHIP_interface.h deploy/

rm -rf ./optimized
rm config.log
