#!/bin/sh

# Compile using all cores
case "$(uname)" in
    Linux)  NCORES=$(grep -c ^processor /proc/cpuinfo) ;;
    Darwin) NCORES=$(sysctl -n hw.ncpu) ;;
    *)      NCORES=4 ;;
esac
echo "compile with $NCORES cores"

for program in \
    parhip \
    graph2binary \
    graph2binary_external \
    toolbox \
    ;
do
    scons program=$program variant=optimized_nooutput -j $NCORES || {
        echo "compile error in $program. exiting."
        exit
    }
done

rm config.log
