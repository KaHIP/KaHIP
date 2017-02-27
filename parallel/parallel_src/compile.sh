#!/bin/bash


for program in parhip graph2binary graph2binary_external toolbox; do 
scons program=$program variant=optimized -j 16
if [ "$?" -ne "0" ]; then 
        echo "compile error in $program. exiting."
        exit
fi
done

rm config.log
