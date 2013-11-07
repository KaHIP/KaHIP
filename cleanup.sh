#!/bin/bash


scons program=kaffpa variant=optimized -j 4 -c 
scons program=kaffpaE variant=optimized -j 4 -c 
scons program=buffoon variant=optimized -j 4 -c
scons program=kabar variant=optimized -j 4 -c
scons program=partition_to_vertex_separator variant=optimized -j 4 -c
scons program=library variant=optimized -j 4 -c

rm -rf deploy
rm -rf optimized
