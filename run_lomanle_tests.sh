#!/bin/bash

# --- 1. Compilation ---
# ./build.sh
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
    echo "Usage: $0 <file.csv|all> <k_min_list> [manifold_dim] [g_threshold] [o_max] [o_min]"
    echo "Example: $0 results/data/arc.csv 10 1 1.0 0.30 0.05"
    exit 1
fi

dataset_input=$1
k_min_list=(${2//,/ })
m_dim=${3:-1}        # Default 1 (curve)
g_thresh=${4:-1.0}   # Default 1.0
o_max=${5:-0.30}     # Default 0.30 (Maximum Overlap)
o_min=${6:-0.05}     # Default 0.05 (Minimum Connectivity)

# --- 3. Processing Function ---
process_file() {
    local f=$1
    for k in "${k_min_list[@]}"; do
        echo "------------------------------------------"
        echo "EXECUTING: $f"
        echo "Parameters: k_min=$k, dim=$m_dim, gap_threshold=$g_thresh, o_max=$o_max, o_min=$o_min"

        # Execute Fortran with the 6 arguments
        ./build/test_lomanle "$f" "$k" "$m_dim" "$g_thresh" "$o_max" "$o_min"

        if [ $? -ne 0 ]; then
            echo "Error during Fortran execution."
            continue
        fi

        # Define output CSV name including k, g, and overlap parameters
        base_name=$(basename "$f" .csv)
        csv_out="results/data/${base_name}_k${k}_g${g_thresh}_omax${o_max}_lomanle.csv"
        mv lomanle_output.csv "$csv_out"

        # R SCRIPT CALL (Passes parameters for plotting)
        echo "Generating visualization..."
        Rscript r/plot_lomanle_spheres.R "$csv_out" "$k" "$g_thresh" "$o_max" "$o_min"
    done
}

# --- 4. Main Logic ---
if [ "$dataset_input" == "all" ]; then
    for f in results/data/*.csv; do
        # Avoid re-processing files that are already lomanle outputs
        if [[ "$f" != *"_lomanle"* ]]; then
            process_file "$f"
        fi
    done
else
    process_file "$dataset_input"
fi

echo "------------------------------------------"
echo "Process finished. Check results/plots/ for PDF reports."