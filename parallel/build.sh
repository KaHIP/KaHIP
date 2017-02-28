#!/bin/bash

rm -rf deploy 2>/dev/null

cd modified_kahip 
./cleanup.sh
./compile.sh
cd ..

cp modified_kahip/deploy/libkahip.a parallel_src/extern/kaHIP_lib/
cp modified_kahip/deploy/kaHIP_interface.h parallel_src/extern/kaHIP_lib/

cd parallel_src/
./cleanup.sh
./compile.sh
cd ..

cp parallel_src/optimized_nooutput/g* ../deploy
cp parallel_src/optimized_nooutput/parhip* ../deploy/parhip
cp parallel_src/optimized_nooutput/toolbox* ../deploy/
