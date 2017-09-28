#!/bin/sh

# Compile using all cores
case "$(uname)" in
    Linux)  NCORES=$(grep -c ^processor /proc/cpuinfo);;
    Darwin) NCORES=$(sysctl -n hw.ncpu) ;;
    *)      NCORES=4 ;;
esac

echo "compile with $NCORES cores"

scons program=library variant=optimized -j $NCORES

rm -rf deploy
mkdir deploy
cp ./optimized/interface/lib* deploy
cp ./interface/kaHIP_interface.h deploy

rm -rf ./optimized
rm config.log
