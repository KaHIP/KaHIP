#!/bin/sh

NCORES=16

for program in \
    parhip \
    edge_list_to_metis_graph \
    friendster_list_to_metis_graph \
    graph2binary \
    graph2binary_external \
    readbgf \
    toolbox \
    ;
do
    scons program=$program variant=optimized -j $NCORES || {
        echo "compile error in $program. exiting."
        exit
    }
done

rm config.log
