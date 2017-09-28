#!/bin/sh
cd ${0%/*} || exit 1   # Run from this directory

rm -rf deploy 2>/dev/null

(
    cd modified_kahip || exit
    ./cleanup.sh
    ./compile.sh
)

mkdir -p parallel_src/extern/kaHIP_lib

for file in \
    libkahip.a \
    kaHIP_interface.h \
    ;
do
    cp modified_kahip/deploy/$file  parallel_src/extern/kaHIP_lib
done

(
    cd parallel_src || exit
    ./cleanup.sh
    ./compile.sh
)

mkdir -p ../deploy

cp parallel_src/optimized_nooutput/g* ../deploy
cp parallel_src/optimized_nooutput/parhip* ../deploy/parhip
cp parallel_src/optimized_nooutput/toolbox* ../deploy

# -----------------------------------------------------------------------------
