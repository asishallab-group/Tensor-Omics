/* stc_cli: a command line interface to Shape Truthful Clustering (STC), see
 * misc/mod_STC.md, "Command line interface (CLI) in C". Reads an all-real-number CSV via
 * csv_table (libcsv-backed, see csv_table.h), runs the STC pipeline through the generated
 * `_c` bindings (never re-implementing any of STC itself -- this file is wiring only), and
 * writes an interactive HTML/D3 report, the same results as JSON, and three CSV/TSV
 * companions for post-processing in Python/R/etc., all into --output-dir.
 *
 * Every STC/estimation parameter that has no compile-time-constant kernel default
 * (`chordal_dist_max_as_prcnt_of_range`, `d_max`, `G_max`, `RMSE_change_max`, `o`) is a
 * required flag. Every parameter that does have one (`k_min`, `k_density`,
 * `bandwidth_percentile`, `exclusion_radius_percentile`, `f_max`, `a`,
 * `min_overlap_coefficient`, `n_anchors`, `seed_max_set_size`, `first_quartile_percentile`)
 * is optional, defaulting to that same value -- the raw `_c` ABI has no notion of "use the
 * Fortran default" for these (see the generated `src/generated/bindings/c/*_c.F90` files:
 * DM_DEFAULT parameters are documented there, but declared non-optional), so this CLI
 * resolves them itself, once, rather than asking Fortran to.
 *
 * --estimate-parameters runs `estimate_stc_parameters` first and applies its own
 * `k_min`/`k_density`/`chordal_dist_max_as_prcnt_of_range`/`G_max`/`d_max` as the actual run
 * parameters (rounding the two that come back real-valued from a median). Supplying any of
 * those five flags together with --estimate-parameters is a validation error: the user picks
 * estimation or manual values, never both. `density_quantile`, estimation's sixth output, has
 * no corresponding run parameter to apply to at all -- it is reported in the JSON/HTML output
 * (`estimated_density_quantile`) and nowhere else.
 */

#include <argp.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "csv_table.h"

/* ==== the STC/`tox_stc_json`/`tox_stc_csv` `_c` symbols this file calls ========================
 * No shared header exists for these anywhere in the codebase -- every consumer (the R/Python
 * bindings included, see src/generated/bindings/{r,c}/*.c) declares its own prototypes for
 * exactly the symbols it needs, matching the exact `_c.F90` signatures under
 * src/generated/bindings/c/. */

extern void build_kd_index_c(const double *points, const int *n_dimensions, const int *n_points,
                             int *kd_indices, const int *dimension_order, int *ierr);

extern void seeds_c(const double *vectors, const int *n_dimensions, const int *n_vectors,
                    const int *kd_indices, const int *dimension_order, const int *k_density,
                    const double *bandwidth_percentile, const double *exclusion_radius_percentile,
                    unsigned char *is_seed_mask, int *ierr);

extern void estimate_stc_parameters_c(const double *vectors, const int *n_dimensions, const int *n_vectors,
                                      const int *kd_indices, const int *dimension_order,
                                      const int *k_density, const double *bandwidth_percentile,
                                      const int *n_anchors, const double *seed_max_set_size,
                                      const double *first_quartile_percentile,
                                      double *estimated_k_min, double *estimated_k_density,
                                      double *estimated_density_quantile,
                                      double *estimated_chordal_dist_max_as_prcnt_of_range,
                                      double *estimated_G_max, double *estimated_d_max, int *ierr);

extern void ensemble_identification_merged_c(
    const double *vectors, const int *n_dimensions, const int *n_vectors, const int *kd_indices,
    const int *dimension_order, const unsigned char *seed_selection_mask, const int *n_selected_seed,
    const int *k_min, const double *chordal_dist_max_as_prcnt_of_range, const int *d_max,
    const double *G_max, const double *RMSE_change_max, const double *f_max, const int *a, const int *o,
    unsigned char *ensemble_masks, int *ensemble_stop_reason, double *ensemble_growth_radii,
    double *ensemble_U_history, double *ensemble_S_history, int *ensemble_d_history,
    double *ensemble_G_history, double *ensemble_mu_history, int *ensemble_k_history,
    unsigned char *ensemble_accepted_history, int *ensemble_member_added_at_step,
    unsigned char *ensemble_low_confidence_masks, double *ensemble_U_first, int *ensemble_d_first,
    int *ierr);

extern void ensemble_reconciliation_c(const unsigned char *ensemble_masks,
                                      const int *ensemble_stop_reason, const int *n_dimensions,
                                      const int *n_vectors, const int *n_ensembles,
                                      const double *ensemble_U_history, const int *ensemble_d_history,
                                      const double *ensemble_S_history, const double *ensemble_mu_history,
                                      const double *ensemble_G_history, const int *ensemble_k_history,
                                      const unsigned char *ensemble_accepted_history, const int *o,
                                      const char *mode,
                                      const double *min_overlap_coefficient,
                                      const unsigned char *report_overlap_coefficient,
                                      const unsigned char *allowed_stop_reasons,
                                      const int *d_min, const int *d_max, const double *var_explained_min,
                                      const int *max_group_size, int *super_ensembles,
                                      int *n_super_ensembles, double *super_ensembles_overlap_coefficient,
                                      unsigned char *eligible, unsigned char *eligible_by_stop_condition,
                                      unsigned char *eligible_by_dimension, unsigned char *eligible_by_var_explained,
                                      int *ierr);

extern void serialize_stc_results_as_json_c(
    const char *filename, const int *filename_strlen, const int *n_dimensions, const int *n_vectors,
    const int *n_selected_seed, const int *o, const int *max_group_size, const int *n_super_ensembles,
    const double *vectors, const char *dim_names, const int *dim_names_strlen,
    const unsigned char *seed_selection_mask, const unsigned char *ensemble_masks,
    const int *ensemble_stop_reason, const double *ensemble_growth_radii, const double *ensemble_U_history,
    const double *ensemble_S_history, const int *ensemble_d_history, const double *ensemble_G_history,
    const double *ensemble_mu_history, const int *ensemble_k_history,
    const unsigned char *ensemble_accepted_history, const int *ensemble_member_added_at_step,
    const unsigned char *ensemble_low_confidence_masks,
    const double *ensemble_U_first, const int *ensemble_d_first,
    const int *super_ensembles, const int *k_min,
    const int *k_density, const double *chordal_dist_max_as_prcnt_of_range, const int *d_max,
    const double *G_max, const double *RMSE_change_max, const double *f_max, const int *a,
    const double *exclusion_radius_percentile, const double *bandwidth_percentile,
    const char *reconciliation_mode, const double *min_overlap_coefficient,
    const unsigned char *allowed_stop_reasons, const int *filter_dim_min, const int *filter_dim_max,
    const double *filter_var_explained_min, const unsigned char *ensemble_eligible,
    const unsigned char *ensemble_eligible_by_stop_condition, const unsigned char *ensemble_eligible_by_dimension,
    const unsigned char *ensemble_eligible_by_var_explained, const int *estimated_k_min,
    const int *estimated_k_density, const double *estimated_density_quantile,
    const double *estimated_chordal_dist_max_as_prcnt_of_range, const double *estimated_G_max,
    const int *estimated_d_max, int *ierr);

extern void write_stc_interactive_html_report_c(
    const char *filename, const int *filename_strlen, const int *n_dimensions, const int *n_vectors,
    const int *n_selected_seed, const int *o, const int *max_group_size, const int *n_super_ensembles,
    const double *vectors, const char *dim_names, const int *dim_names_strlen,
    const unsigned char *seed_selection_mask, const unsigned char *ensemble_masks,
    const int *ensemble_stop_reason, const double *ensemble_growth_radii, const double *ensemble_U_history,
    const double *ensemble_S_history, const int *ensemble_d_history, const double *ensemble_G_history,
    const double *ensemble_mu_history, const int *ensemble_k_history,
    const unsigned char *ensemble_accepted_history, const int *ensemble_member_added_at_step,
    const unsigned char *ensemble_low_confidence_masks,
    const double *ensemble_U_first, const int *ensemble_d_first,
    const int *super_ensembles, const int *k_min,
    const int *k_density, const double *chordal_dist_max_as_prcnt_of_range, const int *d_max,
    const double *G_max, const double *RMSE_change_max, const double *f_max, const int *a,
    const double *exclusion_radius_percentile, const double *bandwidth_percentile,
    const char *reconciliation_mode, const double *min_overlap_coefficient,
    const unsigned char *allowed_stop_reasons, const int *filter_dim_min, const int *filter_dim_max,
    const double *filter_var_explained_min, const unsigned char *ensemble_eligible,
    const unsigned char *ensemble_eligible_by_stop_condition, const unsigned char *ensemble_eligible_by_dimension,
    const unsigned char *ensemble_eligible_by_var_explained, const int *estimated_k_min,
    const int *estimated_k_density, const double *estimated_density_quantile,
    const double *estimated_chordal_dist_max_as_prcnt_of_range, const double *estimated_G_max,
    const int *estimated_d_max, int *ierr);

extern void serialize_stc_points_as_csv_c(const char *filename, const int *filename_strlen,
                                          const int *n_vectors, const int *n_selected_seed,
                                          const int *max_group_size, const int *n_super_ensembles,
                                          const unsigned char *seed_selection_mask,
                                          const unsigned char *ensemble_masks,
                                          const unsigned char *ensemble_low_confidence_masks,
                                          const int *super_ensembles, int *ierr);

extern void serialize_stc_ensemble_overlap_as_csv_c(const char *filename, const int *filename_strlen,
                                                    const int *n_vectors, const int *n_selected_seed,
                                                    const unsigned char *ensemble_masks, int *ierr);

extern void serialize_stc_super_ensembles_as_tsv_c(const char *filename, const int *filename_strlen,
                                                   const int *max_group_size, const int *n_super_ensembles,
                                                   const int *super_ensembles, int *ierr);

/* ==== small helpers ============================================================== */

/* Fills `dst[0..width)` with `src`'s bytes, null-padding the rest -- the layout
 * `tox_conversions::c_char_1d_as_string`/`c_char_2d_as_string` (used by every generated `_c`
 * wrapper that takes a Fortran string) expects: it scans for a null byte within the slot to
 * find the string's real length, falling back to the whole slot if none is found. */
static void fill_padded_string(char *dst, size_t width, const char *src) {
    size_t len = strlen(src);
    if (len > width) len = width;
    memset(dst, 0, width);
    memcpy(dst, src, len);
}

static void die(const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "stc_cli: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    exit(1);
}

static void check_ierr(int ierr, const char *what) {
    if (ierr != 0) die("%s failed (ierr=%d)", what, ierr);
}

/* Comma-separated Stop-Condition names -> allowed_stop_reasons[4] (Fortran-array order: index
 * 0 = STOP_REASON_MAX_SIZE, 1 = REJECTED_AFTER_STABLE, 2 = REJECTED_IMMEDIATELY,
 * 3 = FIXED_POINT -- matching tox_stc_json's own stc_stop_reason_name strings exactly, see
 * misc/mod_STC.md's "Ensemble Reconciliation"). Starts all-allowed, clears an entry for every
 * name that appears; dies on an unrecognized name. `arg` is never longer than a handful of
 * short names in practice, so a fixed stack buffer (not strdup, to avoid a POSIX feature-test
 * macro dependency for something this small) is enough. */
static void parse_allowed_stop_reasons(const char *arg, unsigned char out[4]) {
    static const char *names[4] = {"max_size", "rejected_after_stable", "rejected_immediately",
                                   "fixed_point"};
    char buf[256];
    out[0] = out[1] = out[2] = out[3] = 1;

    if (strlen(arg) >= sizeof(buf)) die("--reconciliation-exclude-stop-reasons: value too long");
    fill_padded_string(buf, sizeof(buf), arg);

    for (char *tok = strtok(buf, ","); tok != NULL; tok = strtok(NULL, ",")) {
        int matched = 0;
        for (int i = 0; i < 4; i++) {
            if (strcmp(tok, names[i]) == 0) {
                out[i] = 0;
                matched = 1;
                break;
            }
        }
        if (!matched) {
            die("--reconciliation-exclude-stop-reasons: unknown Stop Condition '%s' (expected any of "
                "max_size,rejected_after_stable,rejected_immediately,fixed_point)", tok);
        }
    }
}

static char *join_path(const char *dir, const char *name) {
    size_t len = strlen(dir) + 1 + strlen(name) + 1;
    char *out = malloc(len);
    if (out == NULL) die("out of memory building a path");
    snprintf(out, len, "%s/%s", dir, name);
    return out;
}

/* ==== argument parsing ============================================================ */

enum {
    OPT_SEPARATOR = 256,
    OPT_HEADER,
    OPT_N_RECORDS,
    OPT_K_MIN,
    OPT_K_DENSITY,
    OPT_CHORDAL,
    OPT_D_MAX,
    OPT_G_MAX,
    OPT_RMSE_CHANGE_MAX,
    OPT_F_MAX,
    OPT_A,
    OPT_O,
    OPT_EXCLUSION_RADIUS_PERCENTILE,
    OPT_BANDWIDTH_PERCENTILE,
    OPT_RECONCILIATION_MODE,
    OPT_RECONCILIATION_EXCLUDE_STOP_REASONS,
    OPT_FILTER_DIM_MIN,
    OPT_FILTER_DIM_MAX,
    OPT_FILTER_VAR_EXPLAINED_MIN,
    OPT_MIN_OVERLAP_COEFFICIENT,
    OPT_REPORT_OVERLAP_COEFFICIENT,
    OPT_MAX_GROUP_SIZE,
    OPT_ESTIMATE_PARAMETERS,
    OPT_N_ANCHORS,
    OPT_SEED_MAX_SET_SIZE,
    OPT_FIRST_QUARTILE_PERCENTILE,
};

struct arguments {
    /* input */
    const char *input;
    char separator;
    int header;
    int n_records;
    int have_n_records;
    const char *output_dir;

    /* required: no kernel-side default exists to fall back to */
    double chordal_dist_max_as_prcnt_of_range;
    int have_chordal;
    int d_max;
    int have_d_max;
    double G_max;
    int have_G_max;
    double RMSE_change_max;
    int have_RMSE_change_max;
    int o;
    int have_o;

    /* optional: this CLI resolves the same default the kernel documents */
    int k_min;
    int have_k_min;
    int k_density;
    int have_k_density;
    double bandwidth_percentile;
    int have_bandwidth_percentile;
    double exclusion_radius_percentile;
    double f_max;
    int a;
    const char *reconciliation_mode;
    const char *reconciliation_exclude_stop_reasons;
    int have_filter_dim_min;
    int filter_dim_min;
    int have_filter_dim_max;
    int filter_dim_max;
    int have_filter_var_explained_min;
    double filter_var_explained_min;
    double min_overlap_coefficient;
    int report_overlap_coefficient;
    int max_group_size;
    int have_max_group_size;

    /* --estimate-parameters mode */
    int estimate_parameters;
    int n_anchors;
    double seed_max_set_size;
    double first_quartile_percentile;
};

static struct argp_option options[] = {
    {"input", 'i', "FILE", 0, "Input CSV file (required)", 0},
    {"separator", OPT_SEPARATOR, "CHAR", 0, "Column separator (default: ',')", 0},
    {"header", OPT_HEADER, 0, 0, "First line is a header of dimension names (default: off)", 0},
    {"n-records", OPT_N_RECORDS, "N", 0,
     "Number of data rows in the input CSV, excluding the header line if any (required)", 0},
    {"output-dir", 'o', "DIR", 0, "Directory to write report.html/results.json/*.csv/*.tsv into (required)", 0},

    {"k-min", OPT_K_MIN, "N", 0, "Growth radius neighborhood size (default: 30)", 1},
    {"k-density", OPT_K_DENSITY, "N", 0, "Density estimation neighborhood size (default: 30)", 1},
    {"chordal-dist-max-as-prcnt-of-range", OPT_CHORDAL, "FRACTION", 0,
     "Maximum tolerated chordal distance between tangent bases, as a fraction of its own "
     "range (required unless --estimate-parameters)", 1},
    {"d-max", OPT_D_MAX, "N", 0,
     "Maximum tolerated change in intrinsic dimension (required unless --estimate-parameters)", 1},
    {"g-max", OPT_G_MAX, "REAL", 0,
     "Maximum tolerated |log(G_tp1/G_t)| (required unless --estimate-parameters)", 1},
    {"rmse-change-max", OPT_RMSE_CHANGE_MAX, "REAL", 0,
     "Maximum tolerated |log(RMSE_tp1/RMSE_t)| (required)", 1},
    {"f-max", OPT_F_MAX, "FRACTION", 0, "Ensemble size fraction of N above which growth is abandoned "
     "(default: 0.95)", 1},
    {"a", OPT_A, "N", 0, "Minimum accepted-iteration count for a stable rejection (default: 2)", 1},
    {"o", OPT_O, "N", 0, "Trailing observable-history window depth (required)", 1},
    {"exclusion-radius-percentile", OPT_EXCLUSION_RADIUS_PERCENTILE, "PERCENTILE", 0,
     "Seed coverage/exclusion radius percentile (default: 50.0)", 1},
    {"bandwidth-percentile", OPT_BANDWIDTH_PERCENTILE, "PERCENTILE", 0,
     "Density-estimate kernel bandwidth percentile (default: 68.27)", 1},

    {"reconciliation-mode", OPT_RECONCILIATION_MODE, "MODE", 0,
     "One of report|merge_overlap_coefficient|merge_any (default: merge_overlap_coefficient)", 2},
    {"reconciliation-exclude-stop-reasons", OPT_RECONCILIATION_EXCLUDE_STOP_REASONS, "LIST", 0,
     "Comma-separated Stop Conditions to exclude from reconciliation entirely (never merged, never "
     "in overlap_coefficient_matrix -- still reported everywhere else): any of "
     "max_size,rejected_after_stable,rejected_immediately,fixed_point (default: none excluded)", 2},
    {"filter-dim-min", OPT_FILTER_DIM_MIN, "N", 0,
     "Minimum tolerated final intrinsic dimension for reconciliation eligibility, inclusive "
     "(default: no lower bound) -- distinct from --d-max, which bounds dimension *drift* during "
     "growth, not the final dimension's own value", 2},
    {"filter-dim-max", OPT_FILTER_DIM_MAX, "N", 0,
     "Maximum tolerated final intrinsic dimension for reconciliation eligibility, inclusive "
     "(default: no upper bound)", 2},
    {"filter-var-explained-min", OPT_FILTER_VAR_EXPLAINED_MIN, "FRACTION", 0,
     "Minimum tolerated final classical variance explained (sum(tangent eigenvalues) / "
     "(sum(tangent eigenvalues) + normal_error)) for reconciliation eligibility (default: no "
     "filtering)", 2},
    {"min-overlap-coefficient", OPT_MIN_OVERLAP_COEFFICIENT, "FRACTION", 0,
     "Minimum Overlap Coefficient for merge_overlap_coefficient mode (default: 0.9)", 2},
    {"report-overlap-coefficient", OPT_REPORT_OVERLAP_COEFFICIENT, 0, 0,
     "Also compute each super-ensemble's own merge-chain Overlap Coefficient (default: off)", 2},
    {"max-group-size", OPT_MAX_GROUP_SIZE, "N", 0,
     "Maximum ensembles one super-ensemble can hold (default: min(1024, n_ensembles))", 2},

    {"estimate-parameters", OPT_ESTIMATE_PARAMETERS, 0, 0,
     "Estimate k_min/k_density/chordal_dist_max_as_prcnt_of_range/G_max/d_max from the data and use "
     "them as this run's own values, instead of supplying them by hand -- mutually exclusive with "
     "--k-min/--k-density/--chordal-dist-max-as-prcnt-of-range/--g-max/--d-max", 3},
    {"n-anchors", OPT_N_ANCHORS, "N", 0, "Number of estimator anchors (default: 5)", 3},
    {"seed-max-set-size", OPT_SEED_MAX_SET_SIZE, "PERCENTILE", 0,
     "Cap on total estimator-anchor-cloud growth, as a percent of N (default: 5.0)", 3},
    {"first-quartile-percentile", OPT_FIRST_QUARTILE_PERCENTILE, "PERCENTILE", 0,
     "Percentile used to read off the estimated parameters (default: 25.0)", 3},

    {0},
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *args = state->input;

    switch (key) {
        case 'i': args->input = arg; break;
        case 'o': args->output_dir = arg; break;
        case OPT_SEPARATOR:
            args->separator = (strcmp(arg, "\\t") == 0) ? '\t' : arg[0];
            break;
        case OPT_HEADER: args->header = 1; break;
        case OPT_N_RECORDS: args->n_records = atoi(arg); args->have_n_records = 1; break;

        case OPT_K_MIN: args->k_min = atoi(arg); args->have_k_min = 1; break;
        case OPT_K_DENSITY: args->k_density = atoi(arg); args->have_k_density = 1; break;
        case OPT_CHORDAL:
            args->chordal_dist_max_as_prcnt_of_range = atof(arg);
            args->have_chordal = 1;
            break;
        case OPT_D_MAX: args->d_max = atoi(arg); args->have_d_max = 1; break;
        case OPT_G_MAX: args->G_max = atof(arg); args->have_G_max = 1; break;
        case OPT_RMSE_CHANGE_MAX:
            args->RMSE_change_max = atof(arg);
            args->have_RMSE_change_max = 1;
            break;
        case OPT_F_MAX: args->f_max = atof(arg); break;
        case OPT_A: args->a = atoi(arg); break;
        case OPT_O: args->o = atoi(arg); args->have_o = 1; break;
        case OPT_EXCLUSION_RADIUS_PERCENTILE: args->exclusion_radius_percentile = atof(arg); break;
        case OPT_BANDWIDTH_PERCENTILE:
            args->bandwidth_percentile = atof(arg);
            args->have_bandwidth_percentile = 1;
            break;

        case OPT_RECONCILIATION_MODE: args->reconciliation_mode = arg; break;
        case OPT_RECONCILIATION_EXCLUDE_STOP_REASONS: args->reconciliation_exclude_stop_reasons = arg; break;
        case OPT_FILTER_DIM_MIN: args->filter_dim_min = atoi(arg); args->have_filter_dim_min = 1; break;
        case OPT_FILTER_DIM_MAX: args->filter_dim_max = atoi(arg); args->have_filter_dim_max = 1; break;
        case OPT_FILTER_VAR_EXPLAINED_MIN:
            args->filter_var_explained_min = atof(arg);
            args->have_filter_var_explained_min = 1;
            break;
        case OPT_MIN_OVERLAP_COEFFICIENT: args->min_overlap_coefficient = atof(arg); break;
        case OPT_REPORT_OVERLAP_COEFFICIENT: args->report_overlap_coefficient = 1; break;
        case OPT_MAX_GROUP_SIZE: args->max_group_size = atoi(arg); args->have_max_group_size = 1; break;

        case OPT_ESTIMATE_PARAMETERS: args->estimate_parameters = 1; break;
        case OPT_N_ANCHORS: args->n_anchors = atoi(arg); break;
        case OPT_SEED_MAX_SET_SIZE: args->seed_max_set_size = atof(arg); break;
        case OPT_FIRST_QUARTILE_PERCENTILE: args->first_quartile_percentile = atof(arg); break;

        case ARGP_KEY_END:
            if (args->input == NULL) argp_error(state, "--input is required");
            if (!args->have_n_records) argp_error(state, "--n-records is required");
            if (args->output_dir == NULL) argp_error(state, "--output-dir is required");
            if (!args->have_o) argp_error(state, "--o is required");
            if (!args->have_RMSE_change_max) argp_error(state, "--rmse-change-max is required");

            if (args->estimate_parameters) {
                if (args->have_k_min) argp_error(state, "--k-min conflicts with --estimate-parameters");
                if (args->have_k_density) argp_error(state, "--k-density conflicts with --estimate-parameters");
                if (args->have_chordal)
                    argp_error(state, "--chordal-dist-max-as-prcnt-of-range conflicts with --estimate-parameters");
                if (args->have_G_max) argp_error(state, "--g-max conflicts with --estimate-parameters");
                if (args->have_d_max) argp_error(state, "--d-max conflicts with --estimate-parameters");
            } else {
                if (!args->have_chordal)
                    argp_error(state, "--chordal-dist-max-as-prcnt-of-range is required "
                              "(unless --estimate-parameters)");
                if (!args->have_d_max) argp_error(state, "--d-max is required (unless --estimate-parameters)");
                if (!args->have_G_max) argp_error(state, "--g-max is required (unless --estimate-parameters)");
            }
            break;

        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, 0,
                           "Runs Shape Truthful Clustering on an input CSV and writes an interactive "
                           "HTML/D3 report, the same results as JSON, and CSV/TSV companions.",
                           0, 0, 0};

/* ==== the pipeline ================================================================= */

int main(int argc, char **argv) {
    struct arguments args;
    memset(&args, 0, sizeof(args));
    args.separator = ',';
    args.k_min = 30;
    args.k_density = 30;
    args.bandwidth_percentile = 68.27;
    args.exclusion_radius_percentile = 50.0;
    args.f_max = 0.95;
    args.a = 2;
    args.reconciliation_mode = "merge_overlap_coefficient";
    args.min_overlap_coefficient = 0.9;
    args.n_anchors = 5;
    args.seed_max_set_size = 5.0;
    args.first_quartile_percentile = 25.0;

    argp_parse(&argp, argc, argv, 0, 0, &args);

    struct stat sb;
    if (stat(args.output_dir, &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        die("--output-dir '%s' does not exist or is not a directory", args.output_dir);
    }

    /* ---- read the input table -------------------------------------------------- */
    csv_table table;
    if (csv_table_read(args.input, args.separator, args.header, args.n_records, &table) != 0) {
        die("%s", csv_table_last_error());
    }
    int n_dimensions = table.n_dimensions;
    int n_vectors = table.n_records;

    char *dim_names_buf = NULL;
    int dim_names_strlen = 1;
    if (args.header) {
        for (int i = 0; i < n_dimensions; i++) {
            size_t len = strlen(table.dim_names[i]);
            if ((int)len > dim_names_strlen) dim_names_strlen = (int)len;
        }
        dim_names_buf = malloc((size_t)dim_names_strlen * (size_t)n_dimensions);
        if (dim_names_buf == NULL) die("out of memory building dim_names");
        for (int i = 0; i < n_dimensions; i++) {
            fill_padded_string(dim_names_buf + (size_t)i * (size_t)dim_names_strlen,
                               (size_t)dim_names_strlen, table.dim_names[i]);
        }
    } else {
        dim_names_buf = malloc((size_t)dim_names_strlen * (size_t)n_dimensions);
        if (dim_names_buf == NULL) die("out of memory building dim_names");
        for (int i = 0; i < n_dimensions; i++) {
            char name[32];
            snprintf(name, sizeof(name), "dim%d", i + 1);
            if ((int)strlen(name) > dim_names_strlen) {
                /* widen once, from the first oversized generated name -- "dim" + up to 10
                 * digits safely fits within `name`'s own buffer either way */
                dim_names_strlen = (int)strlen(name);
                free(dim_names_buf);
                dim_names_buf = malloc((size_t)dim_names_strlen * (size_t)n_dimensions);
                if (dim_names_buf == NULL) die("out of memory building dim_names");
                for (int j = 0; j < i; j++) {
                    char prev[32];
                    snprintf(prev, sizeof(prev), "dim%d", j + 1);
                    fill_padded_string(dim_names_buf + (size_t)j * (size_t)dim_names_strlen,
                                       (size_t)dim_names_strlen, prev);
                }
            }
            fill_padded_string(dim_names_buf + (size_t)i * (size_t)dim_names_strlen,
                               (size_t)dim_names_strlen, name);
        }
    }

    /* ---- k-d tree, shared by seeds/estimate_stc_parameters/ensemble_identification ---- */
    int *dimension_order = malloc(sizeof(int) * (size_t)n_dimensions);
    int *kd_indices = malloc(sizeof(int) * (size_t)n_vectors);
    if (dimension_order == NULL || kd_indices == NULL) die("out of memory building the k-d tree inputs");
    for (int i = 0; i < n_dimensions; i++) dimension_order[i] = i + 1;

    int ierr = 0;
    build_kd_index_c(table.data, &n_dimensions, &n_vectors, kd_indices, dimension_order, &ierr);
    check_ierr(ierr, "build_kd_index");

    /* ---- optionally estimate this run's own parameters from the data ---------- */
    int have_estimated = 0;
    double estimated_k_min = 0, estimated_k_density = 0, estimated_density_quantile = 0;
    double estimated_chordal = 0, estimated_G_max = 0, estimated_d_max = 0;

    int effective_k_min = args.k_min;
    int effective_k_density = args.k_density;
    double effective_chordal = args.chordal_dist_max_as_prcnt_of_range;
    double effective_G_max = args.G_max;
    int effective_d_max = args.d_max;

    if (args.estimate_parameters) {
        const double *bandwidth_percentile_p =
            args.have_bandwidth_percentile ? &args.bandwidth_percentile : NULL;

        estimate_stc_parameters_c(table.data, &n_dimensions, &n_vectors, kd_indices, dimension_order,
                                  NULL /* k_density: forbidden as a manual input in this mode */,
                                  bandwidth_percentile_p, &args.n_anchors, &args.seed_max_set_size,
                                  &args.first_quartile_percentile, &estimated_k_min, &estimated_k_density,
                                  &estimated_density_quantile, &estimated_chordal, &estimated_G_max,
                                  &estimated_d_max, &ierr);
        check_ierr(ierr, "estimate_stc_parameters");
        have_estimated = 1;

        effective_k_min = (int)lround(estimated_k_min);
        effective_k_density = (int)lround(estimated_k_density);
        effective_chordal = estimated_chordal;
        effective_G_max = estimated_G_max;
        effective_d_max = (int)lround(estimated_d_max);
    }

    /* ---- seeding ---------------------------------------------------------------- */
    unsigned char *is_seed_mask = malloc((size_t)n_vectors);
    if (is_seed_mask == NULL) die("out of memory allocating is_seed_mask");
    seeds_c(table.data, &n_dimensions, &n_vectors, kd_indices, dimension_order, &effective_k_density,
           &args.bandwidth_percentile, &args.exclusion_radius_percentile, is_seed_mask, &ierr);
    check_ierr(ierr, "seeds");

    int n_selected_seed = 0;
    for (int i = 0; i < n_vectors; i++) {
        if (is_seed_mask[i]) n_selected_seed++;
    }

    int max_group_size = args.have_max_group_size ? args.max_group_size
                                                   : (n_selected_seed < 1024 ? n_selected_seed : 1024);
    if (max_group_size < 2) max_group_size = 2;

    /* ---- ensemble identification -------------------------------------------------- */
    unsigned char *ensemble_masks = malloc((size_t)n_vectors * (size_t)n_selected_seed);
    int *ensemble_stop_reason = malloc(sizeof(int) * (size_t)n_selected_seed);
    double *ensemble_growth_radii = malloc(sizeof(double) * (size_t)n_selected_seed);
    double *ensemble_U_history =
        malloc(sizeof(double) * (size_t)n_dimensions * (size_t)n_dimensions * (size_t)args.o * (size_t)n_selected_seed);
    double *ensemble_S_history =
        malloc(sizeof(double) * (size_t)n_dimensions * (size_t)args.o * (size_t)n_selected_seed);
    int *ensemble_d_history = malloc(sizeof(int) * (size_t)args.o * (size_t)n_selected_seed);
    double *ensemble_G_history = malloc(sizeof(double) * (size_t)args.o * (size_t)n_selected_seed);
    double *ensemble_mu_history =
        malloc(sizeof(double) * (size_t)n_dimensions * (size_t)args.o * (size_t)n_selected_seed);
    int *ensemble_k_history = malloc(sizeof(int) * (size_t)args.o * (size_t)n_selected_seed);
    unsigned char *ensemble_accepted_history = malloc((size_t)args.o * (size_t)n_selected_seed);
    int *ensemble_member_added_at_step = malloc(sizeof(int) * (size_t)n_vectors * (size_t)n_selected_seed);
    unsigned char *ensemble_low_confidence_masks = malloc((size_t)n_vectors * (size_t)n_selected_seed);
    double *ensemble_U_first =
        malloc(sizeof(double) * (size_t)n_dimensions * (size_t)n_dimensions * (size_t)n_selected_seed);
    int *ensemble_d_first = malloc(sizeof(int) * (size_t)n_selected_seed);
    if (ensemble_masks == NULL || ensemble_stop_reason == NULL || ensemble_growth_radii == NULL ||
        ensemble_U_history == NULL || ensemble_S_history == NULL || ensemble_d_history == NULL ||
        ensemble_G_history == NULL || ensemble_mu_history == NULL || ensemble_k_history == NULL ||
        ensemble_accepted_history == NULL || ensemble_member_added_at_step == NULL ||
        ensemble_low_confidence_masks == NULL || ensemble_U_first == NULL || ensemble_d_first == NULL) {
        die("out of memory allocating ensemble_identification_merged's outputs");
    }

    ensemble_identification_merged_c(
        table.data, &n_dimensions, &n_vectors, kd_indices, dimension_order, is_seed_mask,
        &n_selected_seed, &effective_k_min, &effective_chordal, &effective_d_max, &effective_G_max,
        &args.RMSE_change_max, &args.f_max, &args.a, &args.o, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_accepted_history,
        ensemble_member_added_at_step, ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first,
        &ierr);
    check_ierr(ierr, "ensemble_identification_merged");

    /* ---- ensemble reconciliation ------------------------------------------------- */
    unsigned char allowed_stop_reasons_buf[4];
    const unsigned char *allowed_stop_reasons_p = NULL;
    if (args.reconciliation_exclude_stop_reasons != NULL) {
        parse_allowed_stop_reasons(args.reconciliation_exclude_stop_reasons, allowed_stop_reasons_buf);
        allowed_stop_reasons_p = allowed_stop_reasons_buf;
    }
    const int *filter_dim_min_p = args.have_filter_dim_min ? &args.filter_dim_min : NULL;
    const int *filter_dim_max_p = args.have_filter_dim_max ? &args.filter_dim_max : NULL;
    const double *filter_var_explained_min_p =
        args.have_filter_var_explained_min ? &args.filter_var_explained_min : NULL;

    int n_super_ensembles_capacity = n_selected_seed * (n_selected_seed - 1);
    if (n_super_ensembles_capacity < 0) n_super_ensembles_capacity = 0;
    int *super_ensembles = malloc(sizeof(int) * (size_t)max_group_size * (size_t)(n_super_ensembles_capacity > 0 ? n_super_ensembles_capacity : 1));
    if (super_ensembles == NULL) die("out of memory allocating super_ensembles");
    int n_super_ensembles = 0;

    /* Reported alongside every ensemble in the JSON regardless of whether reconciliation itself
     * ran at all (it doesn't below n_selected_seed=2) -- default to "eligible", the same no-op
     * every individual filter criterion already falls back to when its own threshold is absent. */
    unsigned char *ensemble_eligible = malloc((size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    unsigned char *ensemble_eligible_by_stop_condition = malloc((size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    unsigned char *ensemble_eligible_by_dimension = malloc((size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    unsigned char *ensemble_eligible_by_var_explained = malloc((size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    if (ensemble_eligible == NULL || ensemble_eligible_by_stop_condition == NULL ||
        ensemble_eligible_by_dimension == NULL || ensemble_eligible_by_var_explained == NULL) {
        die("out of memory allocating reconciliation-eligibility buffers");
    }
    memset(ensemble_eligible, 1, (size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    memset(ensemble_eligible_by_stop_condition, 1, (size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    memset(ensemble_eligible_by_dimension, 1, (size_t)(n_selected_seed > 0 ? n_selected_seed : 1));
    memset(ensemble_eligible_by_var_explained, 1, (size_t)(n_selected_seed > 0 ? n_selected_seed : 1));

    if (n_selected_seed >= 2) {
        double *super_ensembles_overlap_coefficient =
            malloc(sizeof(double) * (size_t)(max_group_size - 1) * (size_t)n_super_ensembles_capacity);
        if (super_ensembles_overlap_coefficient == NULL) die("out of memory allocating overlap coefficient work space");

        char mode_buf[25];
        fill_padded_string(mode_buf, sizeof(mode_buf), args.reconciliation_mode);
        unsigned char report_oc = args.report_overlap_coefficient ? 1 : 0;

        ensemble_reconciliation_c(ensemble_masks, ensemble_stop_reason, &n_dimensions, &n_vectors,
                                  &n_selected_seed, ensemble_U_history, ensemble_d_history, ensemble_S_history,
                                  ensemble_mu_history, ensemble_G_history, ensemble_k_history,
                                  ensemble_accepted_history, &args.o, mode_buf,
                                  &args.min_overlap_coefficient, &report_oc, allowed_stop_reasons_p,
                                  filter_dim_min_p, filter_dim_max_p, filter_var_explained_min_p, &max_group_size,
                                  super_ensembles, &n_super_ensembles, super_ensembles_overlap_coefficient,
                                  ensemble_eligible, ensemble_eligible_by_stop_condition,
                                  ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained,
                                  &ierr);
        check_ierr(ierr, "ensemble_reconciliation");
        free(super_ensembles_overlap_coefficient);
    }

    /* ---- write every output artifact into --output-dir --------------------------- */
    char mode_buf[25];
    fill_padded_string(mode_buf, sizeof(mode_buf), args.reconciliation_mode);

    const int *estimated_k_min_p = NULL, *estimated_k_density_p = NULL, *estimated_d_max_p = NULL;
    const double *estimated_density_quantile_p = NULL, *estimated_chordal_p = NULL, *estimated_G_max_p = NULL;
    int estimated_k_min_i = 0, estimated_k_density_i = 0, estimated_d_max_i = 0;
    if (have_estimated) {
        estimated_k_min_i = (int)lround(estimated_k_min);
        estimated_k_density_i = (int)lround(estimated_k_density);
        estimated_d_max_i = (int)lround(estimated_d_max);
        estimated_k_min_p = &estimated_k_min_i;
        estimated_k_density_p = &estimated_k_density_i;
        estimated_d_max_p = &estimated_d_max_i;
        estimated_density_quantile_p = &estimated_density_quantile;
        estimated_chordal_p = &estimated_chordal;
        estimated_G_max_p = &estimated_G_max;
    }

    char *html_path = join_path(args.output_dir, "report.html");
    int html_path_len = (int)strlen(html_path);
    write_stc_interactive_html_report_c(
        html_path, &html_path_len, &n_dimensions, &n_vectors, &n_selected_seed, &args.o, &max_group_size,
        &n_super_ensembles, table.data, dim_names_buf, &dim_names_strlen, is_seed_mask, ensemble_masks,
        ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history,
        ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first,
        super_ensembles, &effective_k_min, &effective_k_density,
        &effective_chordal, &effective_d_max, &effective_G_max, &args.RMSE_change_max, &args.f_max,
        &args.a, &args.exclusion_radius_percentile, &args.bandwidth_percentile, mode_buf,
        &args.min_overlap_coefficient, allowed_stop_reasons_p, filter_dim_min_p, filter_dim_max_p,
        filter_var_explained_min_p, ensemble_eligible, ensemble_eligible_by_stop_condition,
        ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min_p, estimated_k_density_p,
        estimated_density_quantile_p, estimated_chordal_p, estimated_G_max_p, estimated_d_max_p, &ierr);
    check_ierr(ierr, "write_stc_interactive_html_report");
    free(html_path);

    char *json_path = join_path(args.output_dir, "results.json");
    int json_path_len = (int)strlen(json_path);
    serialize_stc_results_as_json_c(
        json_path, &json_path_len, &n_dimensions, &n_vectors, &n_selected_seed, &args.o, &max_group_size,
        &n_super_ensembles, table.data, dim_names_buf, &dim_names_strlen, is_seed_mask, ensemble_masks,
        ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history,
        ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first,
        super_ensembles, &effective_k_min, &effective_k_density,
        &effective_chordal, &effective_d_max, &effective_G_max, &args.RMSE_change_max, &args.f_max,
        &args.a, &args.exclusion_radius_percentile, &args.bandwidth_percentile, mode_buf,
        &args.min_overlap_coefficient, allowed_stop_reasons_p, filter_dim_min_p, filter_dim_max_p,
        filter_var_explained_min_p, ensemble_eligible, ensemble_eligible_by_stop_condition,
        ensemble_eligible_by_dimension, ensemble_eligible_by_var_explained, estimated_k_min_p, estimated_k_density_p,
        estimated_density_quantile_p, estimated_chordal_p, estimated_G_max_p, estimated_d_max_p, &ierr);
    check_ierr(ierr, "serialize_stc_results_as_json");
    free(json_path);

    char *points_path = join_path(args.output_dir, "points.csv");
    int points_path_len = (int)strlen(points_path);
    serialize_stc_points_as_csv_c(points_path, &points_path_len, &n_vectors, &n_selected_seed,
                                  &max_group_size, &n_super_ensembles, is_seed_mask, ensemble_masks,
                                  ensemble_low_confidence_masks, super_ensembles, &ierr);
    check_ierr(ierr, "serialize_stc_points_as_csv");
    free(points_path);

    char *overlap_path = join_path(args.output_dir, "ensemble_overlap_coefficients.csv");
    int overlap_path_len = (int)strlen(overlap_path);
    serialize_stc_ensemble_overlap_as_csv_c(overlap_path, &overlap_path_len, &n_vectors, &n_selected_seed,
                                            ensemble_masks, &ierr);
    check_ierr(ierr, "serialize_stc_ensemble_overlap_as_csv");
    free(overlap_path);

    char *super_ensembles_path = join_path(args.output_dir, "super_ensembles.tsv");
    int super_ensembles_path_len = (int)strlen(super_ensembles_path);
    serialize_stc_super_ensembles_as_tsv_c(super_ensembles_path, &super_ensembles_path_len, &max_group_size,
                                           &n_super_ensembles, super_ensembles, &ierr);
    check_ierr(ierr, "serialize_stc_super_ensembles_as_tsv");
    free(super_ensembles_path);

    printf("stc_cli: wrote report.html, results.json, points.csv, ensemble_overlap_coefficients.csv, "
           "super_ensembles.tsv to '%s' (%d vectors, %d dimensions, %d ensembles, %d super-ensembles)\n",
           args.output_dir, n_vectors, n_dimensions, n_selected_seed, n_super_ensembles);

    csv_table_free(&table);
    return 0;
}
