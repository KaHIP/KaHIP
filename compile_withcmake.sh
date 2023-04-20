#!/bin/bash
cd "${0%/*}" || exit  # Run from current directory (source directory) or exit

if [[ -z "$NCORES" ]]; then 
    case "$(uname)" in
        Darwin)
            NCORES=$(sysctl -n hw.ncpu)
            ;;
        *)
            NCORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null)
            ;;
    esac
    [ -n "$NCORES" ] || NCORES=4
fi

rm -rf deploy
rm -rf build
mkdir build

if [ "$1" == "BUILDPYTHONMODULE" ]; then
    ADDITIONAL_ARGS="-DBUILDPYTHONMODULE=On"
else 
    ADDITIONAL_ARGS="$1"
fi

(cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release $ADDITIONAL_ARGS && \
    make -j $NCORES)

echo
echo "Copying files into 'deploy'"
mkdir deploy

# Serial
echo
echo "... executables (serial)"
for name in \
    edge_partitioning \
    evaluator \
    global_multisection \
    graphchecker \
    kaffpa \
    label_propagation \
    node_ordering \
    node_separator \
    partition_to_vertex_separator \
;
do
    cp ./build/"$name" deploy/
done
if [[ -f ./build/kaffpaE ]]; then 
    cp ./build/kaffpaE deploy/
fi

# These may or may not exist
for name in \
    fast_node_ordering \
    ilp_improve \
    ilp_exact \
;
do
    if [ -f build/"$name" ]; then
        cp ./build/"$name" deploy/
    fi
done

echo "... libraries (serial)"
mv -f ./build/libkahip_static.a deploy/libkahip.a
cp ./build/libkahip* deploy/

echo "... headers (serial)"
cp ./interface/kaHIP_interface.h deploy/


# Parallel
echo
if [ -d ./build/parallel ]
then
    echo "... executables (parallel)"

    cp ./build/parallel/parallel_src/dsp* ./deploy/distributed_edge_partitioning
    cp ./build/parallel/parallel_src/g* ./deploy
    cp ./build/parallel/parallel_src/parhip* ./deploy/parhip
    cp ./build/parallel/parallel_src/toolbox* ./deploy/

    echo "... libraries (parallel)"
    cp ./build/parallel/parallel_src/libparhip_inter*.a deploy/libparhip.a
    mkdir deploy/parallel
    cp ./build/parallel/modified_kahip/lib*.a deploy/parallel/libkahip.a

    echo "... headers (parallel)"
    cp ./parallel/parallel_src/interface/parhip_interface.h deploy/

else
    echo "Parhip was not built - skipping deployment"
fi

# maybe adapt paths here and python version here
if [ "$1" == "BUILDPYTHONMODULE" ]; then
    cp misc/pymodule/call* deploy/
    cp ./build/kahip.cp* deploy/
fi

echo
echo "Created files in deploy/"
echo =========================
(cd deploy 2>/dev/null && ls -dF *)
echo =========================
echo
echo "Can remove old build directory"
echo

# ------------------------------------------------------------------------
