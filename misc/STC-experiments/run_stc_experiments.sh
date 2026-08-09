#!/bin/bash
# Orchestrates the STC experiments in this directory: builds the project, generates the
# synthetic datasets if they are not there yet, then runs the STC pipeline
# (seeds -> ensemble_identification_merged -> ensemble_reconciliation) and its plots over
# every combination of the given parameter lists. The STC counterpart of branch
# `smoothing`'s run_lomanle_tests.sh -- same "one CSV or 'all', comma-separated parameter
# lists, nested loop" shape, adapted to STC's own parameters and driven by Python bindings
# instead of a compiled Fortran test binary (see run_stc.py's own header for why that is
# possible here and was not there).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# --- 1. Build -------------------------------------------------------------------------
echo "Building tensor-omics..."
(cd "$REPO_ROOT" && ./build.sh) || { echo "Build failed."; exit 1; }

# --- 2. Generate datasets, if missing --------------------------------------------------
if [ ! -d "$SCRIPT_DIR/data" ]; then
    echo "No data/ directory found -- generating datasets..."
    python3 "$SCRIPT_DIR/generate_datasets.py"
fi

# --- 3. Argument parsing ----------------------------------------------------------------
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <file.csv|all> [k_min_list] [alpha_max_deg_list] [d_max_list] [g_max_list] [reconciliation_mode_list] [min_jsi_list] [k_density_list] [exclusion_radius_percentile_list]"
    echo "Example: $0 data/2d/s_curve_2d_noise_low.csv 30 30 1 3.0 merge_jsi 0.1"
    echo "Example with lists: $0 data/2d/s_curve_2d_noise_low.csv 15,30 30,60 1,2 3.0,8.0 merge_jsi,merge_any 0.1"
    echo "Example, all datasets: $0 all"
    echo "Example, k_density and exclusion_radius_percentile swept independently of k_min:"
    echo "  $0 data/2d/s_curve_2d_noise_low.csv 30 30 1 3.0 merge_jsi 0.1 15,30 25,50"
    echo "Note: k_density_list defaults to 30 (not k_min_list) if omitted -- these two are"
    echo "independent SKG-level parameters, see run_stc.py's own --k-density help text."
    exit 1
fi

dataset_input=$1
k_min_list=(${2:-30}); k_min_list=(${k_min_list[@]//,/ })
alpha_max_deg_list=(${3:-30}); alpha_max_deg_list=(${alpha_max_deg_list[@]//,/ })
d_max_list=(${4:-1}); d_max_list=(${d_max_list[@]//,/ })
g_max_list=(${5:-3.0}); g_max_list=(${g_max_list[@]//,/ })
mode_list=(${6:-merge_jsi}); mode_list=(${mode_list[@]//,/ })
min_jsi_list=(${7:-0.1}); min_jsi_list=(${min_jsi_list[@]//,/ })
# Bug fix: these two used to not be forwarded to run_stc.py at all, so every sweep silently
# ran seeds() at its hardcoded default (k_density=30, exclusion_radius_percentile=50.0)
# regardless of what k_min_list/the results filename's "k..." label said -- see
# misc/STC-experiments/README.md.
k_density_list=(${8:-30}); k_density_list=(${k_density_list[@]//,/ })
exclusion_radius_percentile_list=(${9:-50.0}); exclusion_radius_percentile_list=(${exclusion_radius_percentile_list[@]//,/ })

# --- 4. Processing function --------------------------------------------------------------
process_file() {
    local f=$1
    local base_name
    base_name=$(basename "$f" .csv)
    for k in "${k_min_list[@]}"; do
        for alpha in "${alpha_max_deg_list[@]}"; do
            for d in "${d_max_list[@]}"; do
                for g in "${g_max_list[@]}"; do
                    for mode in "${mode_list[@]}"; do
                        for min_jsi in "${min_jsi_list[@]}"; do
                            for kd in "${k_density_list[@]}"; do
                                for erp in "${exclusion_radius_percentile_list[@]}"; do
                                    echo "------------------------------------------"
                                    echo "EXECUTING: $f"
                                    echo "Parameters: k_min=$k, k_density=$kd, alpha_max_deg=$alpha, d_max=$d, g_max=$g, reconciliation_mode=$mode, min_jsi=$min_jsi, exclusion_radius_percentile=$erp"

                                    prefix="$SCRIPT_DIR/results/data/${base_name}_k${k}_kd${kd}_a${alpha}_d${d}_g${g}_${mode}_j${min_jsi}_erp${erp}"
                                    python3 "$SCRIPT_DIR/run_stc.py" "$f" \
                                        --k-min "$k" --k-density "$kd" --alpha-max-deg "$alpha" --d-max "$d" --g-max "$g" \
                                        --reconciliation-mode "$mode" --min-jsi "$min_jsi" \
                                        --exclusion-radius-percentile "$erp" --estimate-parameters \
                                        --out-prefix "$prefix"

                                    Rscript "$SCRIPT_DIR/plot_stc.R" "$prefix"
                                    python3 "$SCRIPT_DIR/render_interactive.py" "$prefix"
                                    mkdir -p "$SCRIPT_DIR/results/plots"
                                    mv "${prefix}.pdf" "$SCRIPT_DIR/results/plots/$(basename "$prefix").pdf"
                                    mv "${prefix}_interactive.html" "$SCRIPT_DIR/results/plots/$(basename "$prefix")_interactive.html"
                                    # Only present for 3+ ambient dimensions, see plot_stc.R's own 3D report section.
                                    if [ -f "${prefix}_3d.pdf" ]; then
                                        mv "${prefix}_3d.pdf" "$SCRIPT_DIR/results/plots/$(basename "$prefix")_3d.pdf"
                                    fi
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
echo "Done. Check $SCRIPT_DIR/results/plots/ for PDF reports."
