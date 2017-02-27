#!/bin/bash

rm -rf deploy

scons program=library variant=optimized -j 4

mkdir deploy
cp ./optimized/interface/lib* deploy/
cp ./interface/kaHIP_interface.h deploy/

rm -rf ./optimized
rm config.log
