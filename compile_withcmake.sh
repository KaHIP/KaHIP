#!/bin/bash

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi
 
rm -rf deploy
rm -rf build
mkdir build
cd build
if [ "$1" == "BUILDPYTHONMODULE" ]; then
cmake ../ -DCMAKE_BUILD_TYPE=Release 
else 
cmake ../ -DCMAKE_BUILD_TYPE=Release $1
fi
make -j $NCORES
cd ..

mkdir deploy
cp ./build/kaffpa deploy/
cp ./build/evaluator deploy/
cp ./build/label_propagation deploy/
cp ./build/kaffpaE deploy/
cp ./build/graphchecker deploy/
cp ./build/partition_to_vertex_separator deploy/
cp ./build/node_separator deploy/
cp ./build/edge_partitioning deploy/
cp ./build/node_ordering deploy/
cp ./build/global_multisection deploy/

if [[ -f "./build/fast_node_ordering" ]]; then 
        cp ./build/fast_node_ordering deploy/ 
fi

if [[ -f "./build/ilp_improve" ]]; then 
        cp ./build/ilp_improve deploy/ 
fi

if [[ -f "./build/ilp_exact" ]]; then 
        cp ./build/ilp_exact deploy/ 
fi


cp ./build/libinterface_static.a deploy/libkahip.a
cp ./build/libinterface.so deploy/libkahip.so
cp ./build/parallel/parallel_src/dsp* ./deploy/distributed_edge_partitioning
cp ./build/parallel/parallel_src/g* ./deploy
cp ./build/parallel/parallel_src/parhip* ./deploy/parhip
cp ./build/parallel/parallel_src/toolbox* ./deploy/
cp ./build/parallel/parallel_src/libparhip_inter*.a deploy/libparhip.a
cp ./interface/kaHIP_interface.h deploy/
cp ./parallel/parallel_src/interface/parhip_interface.h deploy/
cp build/libinter* deploy/

mkdir deploy/parallel
cp ./build/parallel/modified_kahip/lib*.a deploy/parallel/libkahip.a

# maybe adapt paths here and python version here
if [ "$1" == "BUILDPYTHONMODULE" ]; then
cd misc/pymodule/
PYTHONINCLUDEPATH=`python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])"`
g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` kahip.cpp -L../../deploy/ -lkahip -o kahip`python3-config --extension-suffix` -Ipybind11/include -I$PYTHONINCLUDEPATH
cd ..; cd ..;
cp misc/pymodule/call* deploy/
cp misc/pymodule/kahip.cp* deploy/
fi


