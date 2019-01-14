#!/bin/bash


for program in parhip edge_list_to_metis_graph friendster_list_to_metis_graph  graph2binary graph2binary_external readbgf toolbox dspac; do
scons program=$program variant=optimized_nooutput -c
done

rm -rf deploy
rm -rf optimized*
rm config.log
cd interface/
scons variant=optimized -c
rm -rf optimized
cd ..

