#!/bin/sh

# Compile using all cores
case "$(uname)" in
    Linux)  NCORES=$(grep -c ^processor /proc/cpuinfo) ;;
    Darwin) NCORES=$(sysctl -n hw.ncpu) ;;
    *)      NCORES=4 ;;
esac
echo "compile with $NCORES cores"

# -----------------------------------------------------------------------------

for program in \
    node_separator \
    kaffpa \
    evaluator \
    kaffpaE \
    graphchecker \
    label_propagation \
    partition_to_vertex_separator \
    library \
    ;
do
    scons program=$program variant=optimized -j $NCORES || {
        echo "compile error in $program. exiting."
        exit
    }
done

rm -rf deploy
mkdir deploy

# app files
for program in \
    node_separator \
    kaffpa \
    evaluator \
    kaffpaE \
    graphchecker \
    label_propagation \
    partition_to_vertex_separator \
    ;
do
    cp ./optimized/$program deploy
done

# Library files
cp ./optimized/interface/lib* deploy
cp ./interface/kaHIP_interface.h deploy

rm -rf ./optimized
rm config.log

echo "Now building the PARALLEL programs"
( cd parallel && ./build.sh )

# -----------------------------------------------------------------------------
