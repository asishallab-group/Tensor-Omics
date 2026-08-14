#!/bin/bash
# Orchestrates the STC experiments in this directory: builds the project, generates the
# synthetic datasets if they are not there yet, then runs `stc_cli` (see `C-layer/README.md`)
# over every combination of the given parameter lists, one output directory per combination.
# The STC counterpart of branch `smoothing`'s run_lomanle_tests.sh -- same "one CSV or 'all',
# comma-separated parameter lists, nested loop" shape, now driven by the compiled CLI instead
# of `run_stc.py`/`plot_stc.R`/`render_interactive.py` (removed -- see README.md).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# `stc_cli` has no default for these two -- the underlying kernels document none either, see
# `C-layer/README.md` -- so this script fixes them at `run_stc.py`'s own former defaults
# rather than exposing yet more sweep dimensions. Pass them by hand to `stc_cli` directly if
# you need to vary them.
RMSE_CHANGE_MAX=0.4054651081
O=10

# --- 1. Build -------------------------------------------------------------------------
echo "Building tensor-omics..."
(cd "$REPO_ROOT" && ./build.sh) || { echo "Build failed."; exit 1; }

STC_CLI=$(find "$REPO_ROOT/build" -name stc_cli -type f -executable | head -1)
if [ -z "$STC_CLI" ]; then
    echo "Could not find the built stc_cli executable under $REPO_ROOT/build" >&2
    exit 1
fi

# --- 2. Generate datasets, if missing --------------------------------------------------
if [ ! -d "$SCRIPT_DIR/data" ]; then
    echo "No data/ directory found -- generating datasets..."
    python3 "$SCRIPT_DIR/generate_datasets.py"
fi

# --- 3. Argument parsing ----------------------------------------------------------------
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <file.csv|all> [k_min_list] [chordal_dist_max_as_prcnt_of_range_list] [d_max_list] [g_max_list] [reconciliation_mode_list] [min_overlap_coefficient_list] [k_density_list] [exclusion_radius_percentile_list]"
    echo "Example: $0 data/2d/s_curve_2d_noise_low.csv 30 0.5 1 3.0 merge_overlap_coefficient 0.9"
    echo "Example with lists: $0 data/2d/s_curve_2d_noise_low.csv 15,30 0.5,0.7 1,2 3.0,8.0 merge_overlap_coefficient,merge_any 0.9"
    echo "Example, all datasets: $0 all"
    echo "Example, k_density and exclusion_radius_percentile swept independently of k_min:"
    echo "  $0 data/2d/s_curve_2d_noise_low.csv 30 0.5 1 3.0 merge_overlap_coefficient 0.9 15,30 25,50"
    echo "Note: k_density_list defaults to 30 (not k_min_list) if omitted -- these two are"
    echo "independent SKG-level parameters."
    echo "Each combination's five output files (report.html, results.json, points.csv,"
    echo "ensemble_overlap_coefficients.csv, super_ensembles.tsv) land in their own directory"
    echo "under results/, see README.md."
    exit 1
fi

dataset_input=$1
k_min_list=(${2:-30}); k_min_list=(${k_min_list[@]//,/ })
chordal_dist_max_list=(${3:-0.5}); chordal_dist_max_list=(${chordal_dist_max_list[@]//,/ })
d_max_list=(${4:-1}); d_max_list=(${d_max_list[@]//,/ })
g_max_list=(${5:-3.0}); g_max_list=(${g_max_list[@]//,/ })
mode_list=(${6:-merge_overlap_coefficient}); mode_list=(${mode_list[@]//,/ })
min_overlap_coefficient_list=(${7:-0.9}); min_overlap_coefficient_list=(${min_overlap_coefficient_list[@]//,/ })
k_density_list=(${8:-30}); k_density_list=(${k_density_list[@]//,/ })
exclusion_radius_percentile_list=(${9:-50.0}); exclusion_radius_percentile_list=(${exclusion_radius_percentile_list[@]//,/ })

# --- 4. Processing function --------------------------------------------------------------
process_file() {
    local f=$1
    local base_name
    base_name=$(basename "$f" .csv)
    local n_records=$(( $(wc -l < "$f") - 1 ))
    for k in "${k_min_list[@]}"; do
        for chordal in "${chordal_dist_max_list[@]}"; do
            for d in "${d_max_list[@]}"; do
                for g in "${g_max_list[@]}"; do
                    for mode in "${mode_list[@]}"; do
                        for min_overlap_coefficient in "${min_overlap_coefficient_list[@]}"; do
                            for kd in "${k_density_list[@]}"; do
                                for erp in "${exclusion_radius_percentile_list[@]}"; do
                                    echo "------------------------------------------"
                                    echo "EXECUTING: $f"
                                    echo "Parameters: k_min=$k, k_density=$kd, chordal_dist_max_as_prcnt_of_range=$chordal, d_max=$d, g_max=$g, reconciliation_mode=$mode, min_overlap_coefficient=$min_overlap_coefficient, exclusion_radius_percentile=$erp"

                                    out_dir="$SCRIPT_DIR/results/${base_name}_k${k}_kd${kd}_cd${chordal}_d${d}_g${g}_${mode}_oc${min_overlap_coefficient}_erp${erp}"
                                    mkdir -p "$out_dir"
                                    LD_LIBRARY_PATH="$REPO_ROOT/build" "$STC_CLI" \
                                        --input "$f" --header --n-records "$n_records" --output-dir "$out_dir" \
                                        --k-min "$k" --k-density "$kd" \
                                        --chordal-dist-max-as-prcnt-of-range "$chordal" \
                                        --d-max "$d" --g-max "$g" --rmse-change-max "$RMSE_CHANGE_MAX" --o "$O" \
                                        --reconciliation-mode "$mode" --min-overlap-coefficient "$min_overlap_coefficient" \
                                        --exclusion-radius-percentile "$erp"
                                done
                            done
                        done
                    done
                done
            done
        done
    done
}

# --- 5. Main logic -------------------------------------------------------------------------
if [ "$dataset_input" == "all" ]; then
    for f in "$SCRIPT_DIR"/data/2d/*.csv "$SCRIPT_DIR"/data/3d/*.csv; do
        process_file "$f"
    done
else
    process_file "$dataset_input"
fi

echo "------------------------------------------"
echo "Done. Check $SCRIPT_DIR/results/*/report.html for interactive reports."
