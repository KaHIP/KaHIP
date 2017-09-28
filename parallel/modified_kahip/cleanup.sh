#!/bin/sh

NCORES=4
for program in \
    kaffpa \
    kaffpaE \
    partition_to_vertex_separator \
    library \
    graphchecker \
    label_propagation \
    ;
do
    scons program=$program variant=optimized -c -j $NCORES
done

rm -rf deploy
rm -rf optimized
rm config.log
