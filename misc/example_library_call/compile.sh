#!/bin/sh
cd ${0%/*} || exit 1   # Run from this directory

set -x
toplevel=../../
unset recompile

# Files required
for file in \
    libkahip.a \
    kaHIP_interface.h \
    ;
do
    [ -e $toplevel/deploy/$file ] || {
        recompile=true
        break
    }
done

if [ -n "$recompile" ]
then
(
    cd $toplevel && ./compile.sh
)
fi

# Copy files
for file in \
    libkahip.a \
    kaHIP_interface.h \
    ;
do
    cp $toplevel/deploy/$file .
done

scons program=interfacetest variant=optimized -j 8
cp ./optimized/interface_test .
rm -rf ./optimized/
rm config.log
