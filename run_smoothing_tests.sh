#!/bin/bash

# Compile project
./build.sh

# Compile smoothing program
gfortran -O2 -I build -o build/smooth_all test_aux/smooth_all_methods.f90 build/*.o -L/usr/local/lib -llapack -lblas

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <all|dataset_name> <k_neighbors> <n_iters_max>"
    exit 1
fi

k_neighbors=$2
n_iters_max=$3

if [ "$1" == "all" ]; then
    echo "Running smoothing tests on all datasets with k_neighbors=$k_neighbors and n_iters_max=$n_iters_max..."
    # Run test for all generated datasets
    for f in results/data/*.csv; do
        if [[ "$f" != *"_smoothed"* ]]; then
            echo "Processing: $f"
            ./build/smooth_all "$f" "$k_neighbors" "$n_iters_max"
        fi
    done
else
    echo "Running smoothing test on dataset: $1 with k_neighbors=$k_neighbors and n_iters_max=$n_iters_max"
    ./build/smooth_all "$1" "$k_neighbors" "$n_iters_max"
    exit 0
fi

Rscript r/plot_smooth_all.r "$k_neighbors" "$n_iters_max"