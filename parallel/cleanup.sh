#!/bin/sh
cd ${0%/*} || exit 1   # Run from this directory

for dir in parallel_src modified_kahip
do
(
    cd $dir && ./cleanup.sh
)
done

# -----------------------------------------------------------------------------
