#!/bin/sh

NCORES=4

parallel/cleanup.sh

for program in \
    kaffpa \
    kaffpaE \
    partition_to_vertex_separator \
    library variant=optimized \
    graphchecker variant=optimized \
    label_propagation variant=optimized \
    ;
do
    scons program=$program variant=optimized -c -j $NCORES
done

parallel/cleanup.sh

rm -rf deploy
rm -rf optimized
rm -f config.log
rm -f parallel/parallel_src/.sconsign.dblite
rm -f parallel/parallel_src/extern/kaHIP_lib/libkahip.a

# -----------------------------------------------------------------------------
