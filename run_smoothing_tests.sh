#!/bin/bash

# --- 1. Compilation ---
./build.sh
echo "Compiling test_lomanle with LAPACK support..."
gfortran -O2 -I build -o build/test_lomanle \
    test_aux/test_lomanle.f90 \
    build/*.o \
    -Lexternal/lib -lloess_netlib -llapack -lblas

if [ $? -ne 0 ]; then
    echo "Compilation Error. Verify that lomanle_mod.o is in build/."
    exit 1
fi

# --- 2. Argument Validation ---
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <file.csv|all> <k_min_list> [manifold_dim] [g_threshold_list] [o_max_list] [o_min_list] [stability_threshold_list] [scale_factor_list] [max_iterations_list] [relative_conv_tol_list]"
    echo "Example: $0 results/data/arc.csv 10 1 1.0 0.30 0.05 0.90 2.5 50 0.01"
    echo "Example with lists: $0 results/data/arc.csv 10,20,30 1 1.0,2.0 0.30,0.40 0.05,0.10 0.90 2.5 50 0.01"
    exit 1
fi

dataset_input=$1
k_min_list=(${2//,/ })
m_dim=${3:-1}             # Default 1 (curve)
g_thresh_list=(${4:-1.0})
g_thresh_list=(${g_thresh_list//,/ })   # Convert to array
o_max_list=(${5:-0.30})
o_max_list=(${o_max_list//,/ })         # Convert to array
o_min_list=(${6:-0.05})
o_min_list=(${o_min_list//,/ })         # Convert to array
stability_list=(${7:-0.90})             # Tangent-stability threshold, see smoothing_vecinos.md
stability_list=(${stability_list//,/ })
scale_factor_list=(${8:-2.5})           # Local-scale distance-jump multiplier
scale_factor_list=(${scale_factor_list//,/ })
max_iter_list=(${9:-50})                # Convergence loop iteration cap
max_iter_list=(${max_iter_list//,/ })
conv_tol_list=(${10:-0.01})             # Convergence tolerance, as a fraction of median nearest-neighbor distance
conv_tol_list=(${conv_tol_list//,/ })

# --- 3. Processing Function ---
process_file() {
    local f=$1
    for k in "${k_min_list[@]}"; do
        for g_thresh in "${g_thresh_list[@]}"; do
            for o_max in "${o_max_list[@]}"; do
                for o_min in "${o_min_list[@]}"; do
                    for stability in "${stability_list[@]}"; do
                        for scale_factor in "${scale_factor_list[@]}"; do
                            for max_iter in "${max_iter_list[@]}"; do
                                for conv_tol in "${conv_tol_list[@]}"; do
                                    echo "------------------------------------------"
                                    echo "EXECUTING: $f"
                                    echo "Parameters: k_min=$k, dim=$m_dim, gap_threshold=$g_thresh, o_max=$o_max, o_min=$o_min, stability=$stability, scale_factor=$scale_factor, max_iter=$max_iter, conv_tol=$conv_tol"

                                    # Execute Fortran with the 10 arguments
                                    ./build/test_lomanle "$f" "$k" "$m_dim" "$g_thresh" "$o_max" "$o_min" "$stability" "$scale_factor" "$max_iter" "$conv_tol"

                                    if [ $? -ne 0 ]; then
                                        echo "Error during Fortran execution."
                                        continue
                                    fi

                                    # Define output CSV name including k, g, overlap, stability and iteration parameters
                                    base_name=$(basename "$f" .csv)
                                    csv_out="results/data/${base_name}_k${k}_g${g_thresh}_omax${o_max}_st${stability}_sf${scale_factor}_mi${max_iter}_ct${conv_tol}_lomanle.csv"
                                    edges_out="${csv_out%.csv}_edges.csv"
                                    echo "Saving output to: $csv_out"
                                    mv lomanle_output.csv "$csv_out"
                                    echo "Saving skeleton edges to: $edges_out"
                                    mv lomanle_edges.csv "$edges_out"

                                    # R SCRIPT CALL (Passes parameters for plotting)
                                    echo "Generating visualization..."
                                    Rscript r/plot_lomanle_spheres.R "$csv_out" "$k" "$g_thresh" "$o_max" "$o_min" "$stability" "$scale_factor" "$max_iter" "$conv_tol"
                                done
                            done
                        done
                    done
                done
            done
        done
    done
}

# --- 4. Main Logic ---
if [ "$dataset_input" == "all" ]; then
    for f in results/data/3d/*.csv; do
        # Avoid re-processing files that are already lomanle outputs
        if [[ "$f" != *"_lomanle"* ]]; then
            if [[ "$f" != *"_smoothed"* ]]; then
                process_file "$f"
            fi
        fi
    done
else
    process_file "$dataset_input"
fi

echo "------------------------------------------"
echo "Process finished. Check results/plots/ for PDF reports."