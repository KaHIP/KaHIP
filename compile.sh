#!/bin/bash

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

rm -rf deploy
for program in node_separator kaffpa evaluator kaffpaE graphchecker label_propagation partition_to_vertex_separator library ; do 
scons program=$program variant=optimized -j $NCORES 
if [ "$?" -ne "0" ]; then 
        echo "compile error in $program. exiting."
        exit
fi
done

mkdir deploy
cp ./optimized/kaffpa deploy/
cp ./optimized/evaluator deploy/
cp ./optimized/label_propagation deploy/
cp ./optimized/kaffpaE deploy/
cp ./optimized/graphchecker deploy/
cp ./optimized/partition_to_vertex_separator deploy/
cp ./optimized/interface/lib* deploy/
cp ./optimized/node_separator deploy/
cp ./interface/kaHIP_interface.h deploy/

rm -rf ./optimized
rm config.log


echo "Now building the PARALLEL programs"
cd parallel
./build.sh
cd ..
