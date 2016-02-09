#!/bin/bash

rm -rf deploy

scons program=kaffpa variant=optimized -j 4 
scons program=evaluator variant=optimized -j 4 
scons program=kaffpaE variant=optimized -j 4
scons program=graphchecker variant=optimized -j 4
scons program=label_propagation variant=optimized -j 4
scons program=partition_to_vertex_separator variant=optimized -j 4
scons program=library variant=optimized -j 4

mkdir deploy
cp ./optimized/kaffpa deploy/
cp ./optimized/evaluator deploy/
cp ./optimized/label_propagation deploy/
cp ./optimized/kaffpaE deploy/
cp ./optimized/graphchecker deploy/
cp ./optimized/partition_to_vertex_separator deploy/
cp ./optimized/interface/lib* deploy/
cp ./interface/kaHIP_interface.h deploy/

rm -rf ./optimized
