#!/bin/bash

rm -rf deploy
for program in node_separator kaffpa evaluator kaffpaE graphchecker label_propagation partition_to_vertex_separator library ; do 
scons program=$program variant=optimized -j 4 
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
