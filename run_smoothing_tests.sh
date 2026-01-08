#!/bin/bash

# Compile project
./build.sh

# Compile smoothing program
gfortran -O2 -I build -o build/smooth_all test_aux/smooth_all_methods.f90 build/*.o -L/usr/local/lib -llapack -lblas

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <all|dataset_name> <k_neighbors_list> <n_iters_max_list>"
    exit 1
fi

k_neighbors_list=(${2//,/ })  # Split comma-separated list into array
n_iters_max_list=(${3//,/ })  # Split comma-separated list into array
method_id=$4

if [ "$1" == "all" ]; then
    echo "Running smoothing tests on all datasets with method_id=$method_id..."
    # Run test for all generated datasets
    for f in results/data/*.csv; do
        if [[ "$f" != *"_smoothed"* ]]; then
            echo "Processing: $f"
            for k_neighbors in "${k_neighbors_list[@]}"; do
                for n_iters_max in "${n_iters_max_list[@]}"; do
                    echo "Running with k_neighbors=$k_neighbors, n_iters_max=$n_iters_max"
                    ./build/smooth_all "$f" "$k_neighbors" "$n_iters_max" "$method_id"
                done
            done
        fi
    done
else
    echo "Running smoothing test on dataset: $1 with method_id=$method_id..."
    for k_neighbors in "${k_neighbors_list[@]}"; do
        for n_iters_max in "${n_iters_max_list[@]}"; do
            echo "Running with k_neighbors=$k_neighbors, n_iters_max=$n_iters_max"
            ./build/smooth_all "$1" "$k_neighbors" "$n_iters_max" "$method_id"
        done
    done
    exit 0
fi

for k_neighbors in "${k_neighbors_list[@]}"; do
    for n_iters_max in "${n_iters_max_list[@]}"; do
        Rscript r/plot_smooth_all.r "$k_neighbors" "$n_iters_max"
    done
done

for k_neighbors in "${k_neighbors_list[@]}"; do
    for n_iters_max in "${n_iters_max_list[@]}"; do
        Rscript r/plot_anwil_std.R "$k_neighbors" "$n_iters_max"
    done
done