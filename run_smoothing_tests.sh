#!/bin/bash

# Compile project
# ./build.sh

# # Compile smoothing program
# gfortran -O2 -I build -o build/smooth_all test_aux/smooth_all_methods.f90 build/*.o -Lexternal/lib -lloess_netlib -llapack -lblas

# if [ "$#" -lt 3 ]; then
#     echo "Usage: $0 <all|dataset_name> <k_neighbors_list> <span_list> <n_iters_max_list> <kernel_type_list> <k_neighbors_sigma_list> <method_id>"
#     exit 1
# fi

k_neighbors_list=(${2//,/ })  # Split comma-separated list into array
span_list=(${3//,/ })  # Split comma-separated list into array
n_iters_max_list=(${4//,/ })  # Split comma-separated list into array
kernel_type_list=(${5//,/ })  # Split comma-separated list into array
k_neighbors_sigma_list=(${6//,/ })  # Split comma-separated list into array
method_id=$7

# if [ "$1" == "all" ]; then
#     echo "Running smoothing tests on all datasets with method_id=$method_id..."
#     # Run test for all generated datasets
#     for f in results/data/*.csv; do
#         if [[ "$f" != *"_smoothed"* ]]; then
#             echo "Processing: $f"
#             for k_neighbors in "${k_neighbors_list[@]}"; do
#                 for span in "${span_list[@]}"; do
#                     for n_iters_max in "${n_iters_max_list[@]}"; do
#                         for kernel_type in "${kernel_type_list[@]}"; do
#                             for k_neighbors_sigma in "${k_neighbors_sigma_list[@]}"; do
#                                 echo "Running with k_neighbors=$k_neighbors, span=$span, n_iters_max=$n_iters_max, kernel_type=$kernel_type, k_neighbors_sigma=$k_neighbors_sigma"
#                                 ./build/smooth_all "$f" "$k_neighbors" "$n_iters_max" "$method_id" "$k_neighbors_sigma" "$kernel_type" "$span" 
#                             done
#                         done
#                     done
#                 done
#             done
#         fi
#     done
# else
#     echo "Running smoothing test on dataset: $1 with method_id=$method_id..."
#     for k_neighbors in "${k_neighbors_list[@]}"; do
#         for span in "${span_list[@]}"; do
#             for n_iters_max in "${n_iters_max_list[@]}"; do
#                 for kernel_type in "${kernel_type_list[@]}"; do
#                     for k_neighbors_sigma in "${k_neighbors_sigma_list[@]}"; do
#                         echo "Running with k_neighbors=$k_neighbors, span=$span, n_iters_max=$n_iters_max, kernel_type=$kernel_type, k_neighbors_sigma=$k_neighbors_sigma"
#                         ./build/smooth_all "$1" "$k_neighbors" "$n_iters_max" "$method_id" "$k_neighbors_sigma" "$kernel_type" "$span" 
#                     done
#                 done
#             done
#         done
#     done
# fi

for k_neighbors in "${k_neighbors_list[@]}"; do
    echo "Entering k_neighbors loop with k_neighbors=$k_neighbors"
    for span in "${span_list[@]}"; do
        echo "  Entering span loop with span=$span"
        for n_iters_max in "${n_iters_max_list[@]}"; do
            for kernel_type in "${kernel_type_list[@]}"; do
                for k_neighbors_sigma in "${k_neighbors_sigma_list[@]}"; do
                    echo "    Entering n_iters_max loop with n_iters_max=$n_iters_max, kernel_type=$kernel_type, k_neighbors_sigma=$k_neighbors_sigma"
                    Rscript r/plot_smooth_all.r "$k_neighbors" "$span" "$n_iters_max" "$kernel_type" "$k_neighbors_sigma"
                done
            done
        done
    done
done

for k_neighbors in "${k_neighbors_list[@]}"; do
    for span in "${span_list[@]}"; do
        for n_iters_max in "${n_iters_max_list[@]}"; do
            for kernel_type in "${kernel_type_list[@]}"; do
                for k_neighbors_sigma in "${k_neighbors_sigma_list[@]}"; do
                    Rscript r/plot_anwil_std.R "$k_neighbors" "$span" "$n_iters_max" "$kernel_type" "$k_neighbors_sigma"
                done
            done
        done
    done
done
