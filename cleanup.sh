#!/bin/bash

cd parallel/parallel_src/
./cleanup.sh
cd ../..
cd parallel/modified_kahip/
./cleanup.sh
cd ../../

scons program=kaffpa variant=optimized -j 4 -c 
scons program=kaffpaE variant=optimized -j 4 -c 
scons program=partition_to_vertex_separator variant=optimized -j 4 -c
scons program=library variant=optimized -j 4 -c
scons program=graphchecker variant=optimized -j 4 -c
scons program=label_propagation variant=optimized -j 4 -c

cd parallel/parallel_src/
./cleanup.sh
cd ../../
cd parallel/modified_kahip/
./cleanup.sh
cd ../../
rm -rf deploy
rm -rf optimized
rm config.log
rm parallel/parallel_src/.sconsign.dblite
rm parallel/parallel_src/extern/kaHIP_lib/libkahip.a

