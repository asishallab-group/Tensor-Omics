# Generated. Do not edit.

#' Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
#'
#' The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
#'
#' 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
#' the caller, e.g. with
#' \code{\link{pool_means}}) form its
#' neighborhoods
#' (\code{\link{construct_neighborhoods_ranged}}).
#' 2. Per reference point, all studies' residuals are pooled and the bin-count search
#' `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
#' \code{\link{determine_bin_count_occupancy}}
#' by default, or Stone's rule
#' \code{\link{determine_bin_count_stone}}).
#' 3. Each study's histogram and pmf is built per point
#' (\code{\link{build_residual_histograms}}).
#' 4. The consensus pmf
#' (\code{\link{create_mean_pmf}}), each
#' study's per-point JSD against it and its weighted global JSD, and the permutation test
#' (\code{\link{gjct_permutation_test}}) give
#' `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
#' `n_permutations = 0` both p-values are NaN.
#'
#' There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
#' still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
#' caller that needs such neighborhoods rejected should use
#' \code{\link{run_js_comp_test_parameter_search}}.
#'
#' The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
#' second extent \code{MAX_N_BINS}
#' (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
#' zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
#' `[, seq_len(max_n_bins_per_point), ]`) itself.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
#' reference point (for example a study whose gene means are all missing), since the
#' equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
#' then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
#' generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
#' values from before the permutation test. An out-of-range neighbor index is also reported
#' as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
#'
#' The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
#' last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
#' caller instead. That order must still be ascending with NaN last (only the order among tied
#' genes is the caller's choice); it is checked neither to be a permutation of
#' `1..max_n_genes_all_studies` per column nor to be sorted.
#' The allocating wrapper runs \code{prepare_gene_means_perms} first, which may handle the call and skip this one.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_run::run_js_comp_test}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point computes \code{gene_means_perms}, \code{tmp_neighborhood_residuals_gathered}, \code{tmp_permutation_mean_pmf}, \code{tmp_permutation_point_pool}, \code{tmp_permutation_point_draw}, \code{tmp_pooled_residuals}, \code{tmp_pooled_residuals_perm}, \code{tmp_bin_counts_search}, \code{tmp_residual_to_study_perm}, \code{tmp_study_block_ends}, \code{tmp_study_block_sizes}, \code{tmp_point_n_neighbors}, \code{tmp_permutation_mean_pmf_counts}, \code{tmp_permutation_counts}, \code{tmp_permutation_pmfs}, \code{tmp_permutation_js_divergences}, \code{tmp_permutation_weights} and \code{tmp_permutation_global_js_divergence} for you and runs \code{prepare_gene_means_perms} first, which may answer the call outright.
#' Call \code{run_js_comp_test_expert} to do that yourself.
#'
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param gene_means a numeric matrix. Per-gene mean expression values for all studies
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study
#'   NaN is permitted for this value.
#' @param x_star a numeric vector. Mean-expression reference points
#' @param n_permutations a integer scalar. Number of permutations of the permutation test; `0` skips it (p-values NaN)
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a bin count to be
#'   admissible, in every `bin_count_mode`
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer array of rank 3. Gene indices of the selected neighborhood, per reference point, per study}
#'   \item{neighborhood_range}{a integer array of rank 3. For each reference point and study, the `[min_idx, max_idx]` neighborhood span}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's histogram bin count `M_j`, from the bin-count search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the bin-count search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the bin-count search}
#'   \item{max_n_bins_per_point}{a integer scalar. `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
#'     `mean_pmf` and `mean_pmf_counts` (see above)}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff the bin-count search failed for this point; it then uses `m_min` bins and
#'     still contributes (see above)}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's non-NaN pooled residual count (N_j)}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic (not used for the decision)}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
#'     are meaningful}
#'   \item{counts}{a integer array of rank 3. Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN replicates per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
#'     meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates per reference point of the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, the primary result; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
#'     comparison only; NaN when `n_permutations` is 0}
#' @export
run_js_comp_test <- function(n_neighbors, gene_means, residuals, x_star, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, bin_count_mode = "bin_count_occupancy") {
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    x_star <- .tox_as_double_vector(x_star, "x_star")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_call", n_neighbors, gene_means, residuals, x_star, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_indices = .result$neighborhood_indices,
        neighborhood_range = .result$neighborhood_range,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        max_n_bins_per_point = .result$max_n_bins_per_point,
        occupancy_failed = .result$occupancy_failed,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins,
        pmfs = .result$pmfs,
        counts = .result$counts,
        included_n_reps = .result$included_n_reps,
        mean_pmf = .result$mean_pmf,
        mean_pmf_counts = .result$mean_pmf_counts,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        js_divergences = .result$js_divergences,
        weights = .result$weights,
        global_js_divergence = .result$global_js_divergence,
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
#'
#' The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
#'
#' 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
#' the caller, e.g. with
#' \code{\link{pool_means}}) form its
#' neighborhoods
#' (\code{\link{construct_neighborhoods_ranged}}).
#' 2. Per reference point, all studies' residuals are pooled and the bin-count search
#' `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
#' \code{\link{determine_bin_count_occupancy}}
#' by default, or Stone's rule
#' \code{\link{determine_bin_count_stone}}).
#' 3. Each study's histogram and pmf is built per point
#' (\code{\link{build_residual_histograms}}).
#' 4. The consensus pmf
#' (\code{\link{create_mean_pmf}}), each
#' study's per-point JSD against it and its weighted global JSD, and the permutation test
#' (\code{\link{gjct_permutation_test}}) give
#' `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
#' `n_permutations = 0` both p-values are NaN.
#'
#' There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
#' still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
#' caller that needs such neighborhoods rejected should use
#' \code{\link{run_js_comp_test_parameter_search}}.
#'
#' The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
#' second extent \code{MAX_N_BINS}
#' (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
#' zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
#' `[, seq_len(max_n_bins_per_point), ]`) itself.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
#' reference point (for example a study whose gene means are all missing), since the
#' equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
#' then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
#' generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
#' values from before the permutation test. An out-of-range neighbor index is also reported
#' as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
#'
#' The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
#' last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
#' caller instead. That order must still be ascending with NaN last (only the order among tied
#' genes is the caller's choice); it is checked neither to be a permutation of
#' `1..max_n_genes_all_studies` per column nor to be sorted.
#' The allocating wrapper runs \code{prepare_gene_means_perms} first, which may handle the call and skip this one.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_run::run_js_comp_test_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{gene_means_perms}, \code{tmp_neighborhood_residuals_gathered}, \code{tmp_permutation_mean_pmf}, \code{tmp_permutation_point_pool}, \code{tmp_permutation_point_draw}, \code{tmp_pooled_residuals}, \code{tmp_pooled_residuals_perm}, \code{tmp_bin_counts_search}, \code{tmp_residual_to_study_perm}, \code{tmp_study_block_ends}, \code{tmp_study_block_sizes}, \code{tmp_point_n_neighbors}, \code{tmp_permutation_mean_pmf_counts}, \code{tmp_permutation_counts}, \code{tmp_permutation_pmfs}, \code{tmp_permutation_js_divergences}, \code{tmp_permutation_weights} and \code{tmp_permutation_global_js_divergence} yourself.
#' \code{run_js_comp_test} computes \code{gene_means_perms}, \code{tmp_neighborhood_residuals_gathered}, \code{tmp_permutation_mean_pmf}, \code{tmp_permutation_point_pool}, \code{tmp_permutation_point_draw}, \code{tmp_pooled_residuals}, \code{tmp_pooled_residuals_perm}, \code{tmp_bin_counts_search}, \code{tmp_residual_to_study_perm}, \code{tmp_study_block_ends}, \code{tmp_study_block_sizes}, \code{tmp_point_n_neighbors}, \code{tmp_permutation_mean_pmf_counts}, \code{tmp_permutation_counts}, \code{tmp_permutation_pmfs}, \code{tmp_permutation_js_divergences}, \code{tmp_permutation_weights} and \code{tmp_permutation_global_js_divergence} for you and runs \code{prepare_gene_means_perms} first, which may answer the call outright.
#'
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param gene_means a numeric matrix. Per-gene mean expression values for all studies
#'   NaN is permitted for this value.
#' @param gene_means_perms a integer matrix. Per study (column), a permutation of that study's `gene_means` in ascending order, NaN last.
#'   The plain entry point builds it for you (heapsort); the expert one takes yours, which
#'   must be sorted that way too (only the order among tied genes is free). It is checked
#'   neither to be a permutation of `1..max_n_genes_all_studies` per column nor to be sorted.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study
#'   NaN is permitted for this value.
#' @param x_star a numeric vector. Mean-expression reference points
#' @param n_permutations a integer scalar. Number of permutations of the permutation test; `0` skips it (p-values NaN)
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a bin count to be
#'   admissible, in every `bin_count_mode`
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{neighborhood_indices}{a integer array of rank 3. Gene indices of the selected neighborhood, per reference point, per study}
#'   \item{neighborhood_range}{a integer array of rank 3. For each reference point and study, the `[min_idx, max_idx]` neighborhood span}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's histogram bin count `M_j`, from the bin-count search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the bin-count search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the bin-count search}
#'   \item{max_n_bins_per_point}{a integer scalar. `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
#'     `mean_pmf` and `mean_pmf_counts` (see above)}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff the bin-count search failed for this point; it then uses `m_min` bins and
#'     still contributes (see above)}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's non-NaN pooled residual count (N_j)}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic (not used for the decision)}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
#'     are meaningful}
#'   \item{counts}{a integer array of rank 3. Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN replicates per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
#'     meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN replicates per reference point of the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, the primary result; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
#'     comparison only; NaN when `n_permutations` is 0}
#' @export
run_js_comp_test_expert <- function(n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, bin_count_mode = "bin_count_occupancy") {
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    gene_means_perms <- .tox_as_integer_matrix(gene_means_perms, "gene_means_perms")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    x_star <- .tox_as_double_vector(x_star, "x_star")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(gene_means_perms)[2] != dim(gene_means)[2])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[2], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(gene_means_perms)[1] != dim(gene_means)[1])
        .tox_shape_error("gene_means_perms", dim(gene_means_perms)[1], "gene_means", dim(gene_means)[1])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_expert_call", n_neighbors, gene_means, gene_means_perms, residuals, x_star, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "tmp_neighborhood_residuals_gathered", "tmp_permutation_mean_pmf", "tmp_permutation_point_pool", "tmp_permutation_point_draw", "tmp_pooled_residuals", "tmp_pooled_residuals_perm", "tmp_bin_counts_search", "tmp_residual_to_study_perm", "tmp_study_block_ends", "tmp_study_block_sizes", "tmp_point_n_neighbors", "tmp_permutation_mean_pmf_counts", "tmp_permutation_counts", "tmp_permutation_pmfs", "tmp_permutation_js_divergences", "tmp_permutation_weights", "tmp_permutation_global_js_divergence", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighborhood_indices = .result$neighborhood_indices,
        neighborhood_range = .result$neighborhood_range,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        max_n_bins_per_point = .result$max_n_bins_per_point,
        occupancy_failed = .result$occupancy_failed,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins,
        pmfs = .result$pmfs,
        counts = .result$counts,
        included_n_reps = .result$included_n_reps,
        mean_pmf = .result$mean_pmf,
        mean_pmf_counts = .result$mean_pmf_counts,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        js_divergences = .result$js_divergences,
        weights = .result$weights,
        global_js_divergence = .result$global_js_divergence,
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}

#' Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test}}: the
#' same pipeline and outputs, on neighborhoods whose size varies per reference point and per
#' study. Call
#' \code{\link{construct_adaptive_neighborhoods}}
#' first and pass its `pooled_neighborhood_range` (its first `n_points` columns) with the same
#' `gene_means` and `residuals`. Each range is decoded back to per-study gene lists; their sizes
#' are returned as `n_neighbors_per_point`. Neighborhoods that decode to the same gene sets as a
#' fixed-k run give bit-identical results. There is no admissibility gate.
#'
#' Errors (`ierr`): `ERR_INVALID_INPUT` before any histogram work for the first range that is
#' reversed, reaches past the non-NaN pooled means, or leaves a study without a gene (every
#' output then undefined); otherwise as in `run_js_comp_test_impl` (a study without a non-NaN
#' residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).
#'
#' The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
#' values each, the size of `residuals`, since one study may contribute all its genes to one
#' neighborhood. The bin-sized arrays have the fixed second extent
#' \code{MAX_N_BINS}; only their first
#' `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
#' `[, seq_len(max_n_bins_per_point), ]`).
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_run::run_js_comp_test_adaptive}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene; the array the
#'   neighborhoods were constructed from
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param pooled_neighborhood_range a integer matrix. For each neighborhood, its first and last sorted pooled position, from
#'   construct_adaptive_neighborhoods
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies*n_studies`.
#' @param n_permutations a integer scalar. Number of permutations of the permutation test; `0` skips it (p-values NaN)
#'   The minimum valid value is `0`.
#'   The default value is `1000`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach for a bin count to be
#'   admissible, in every `bin_count_mode`
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{n_neighbors_per_point}{a integer matrix. For each neighborhood, how many of its genes belong to each study (at least 1 each)}
#'   \item{n_bins_per_point}{a integer vector. Each reference point's selected histogram bin count, from the bin-count search}
#'   \item{shared_residual_range_low}{a numeric vector. Each reference point's lower residual-range bound (R_low), from the bin-count search}
#'   \item{shared_residual_range_high}{a numeric vector. Each reference point's upper residual-range bound (R_high), from the bin-count search}
#'   \item{max_n_bins_per_point}{a integer scalar. `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
#'     `mean_pmf` and `mean_pmf_counts`}
#'   \item{occupancy_failed}{a logical vector. `TRUE` iff the bin-count search failed for this point; it then uses `m_min` bins and
#'     still contributes}
#'   \item{n_pooled_residuals}{a integer vector. Each reference point's non-NaN pooled residual count (N_j)}
#'   \item{min_bin_occupancy}{a integer vector. Each reference point's minimum bin occupancy at `n_bins_per_point`}
#'   \item{mean_bin_occupancy}{a numeric vector. Each reference point's mean bin occupancy at `n_bins_per_point`}
#'   \item{max_bin_occupancy}{a integer vector. Each reference point's maximum bin occupancy at `n_bins_per_point`}
#'   \item{sturges_bins}{a integer vector. Each reference point's Sturges' rule diagnostic (not used for the decision)}
#'   \item{fd_bins}{a integer vector. Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)}
#'   \item{pmfs}{a numeric array of rank 3. `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
#'     are meaningful}
#'   \item{counts}{a integer array of rank 3. Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{included_n_reps}{a integer matrix. Count of non-NaN residuals binned per reference point, per study}
#'   \item{mean_pmf}{a numeric matrix. The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful}
#'   \item{mean_pmf_counts}{a integer matrix. Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
#'     meaningful}
#'   \item{mean_pmf_included_n_reps}{a integer vector. Count of non-NaN residuals per reference point for the consensus pmf}
#'   \item{js_divergences}{a numeric matrix. Per-reference-point JSD of each study against the consensus pmf}
#'   \item{weights}{a numeric matrix. Per-reference-point weights for `global_js_divergence`}
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus pmf}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, the primary result; NaN when `n_permutations` is 0}
#'   \item{p_values_observed_consensus}{a numeric vector. Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
#'     comparison only; NaN when `n_permutations` is 0}
#' @export
run_js_comp_test_adaptive <- function(gene_means, residuals, pooled_neighborhood_range, n_permutations = 1000L, random_seed = 42L, min_residuals_per_bin = 10L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    pooled_neighborhood_range <- .tox_as_integer_matrix(pooled_neighborhood_range, "pooled_neighborhood_range")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_call", gene_means, residuals, pooled_neighborhood_range, n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "gene_means", "residuals", "pooled_neighborhood_range", "n_neighbors_per_point", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", "pooled_neighborhood_range", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_neighbors_per_point = .result$n_neighbors_per_point,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        max_n_bins_per_point = .result$max_n_bins_per_point,
        occupancy_failed = .result$occupancy_failed,
        n_pooled_residuals = .result$n_pooled_residuals,
        min_bin_occupancy = .result$min_bin_occupancy,
        mean_bin_occupancy = .result$mean_bin_occupancy,
        max_bin_occupancy = .result$max_bin_occupancy,
        sturges_bins = .result$sturges_bins,
        fd_bins = .result$fd_bins,
        pmfs = .result$pmfs,
        counts = .result$counts,
        included_n_reps = .result$included_n_reps,
        mean_pmf = .result$mean_pmf,
        mean_pmf_counts = .result$mean_pmf_counts,
        mean_pmf_included_n_reps = .result$mean_pmf_included_n_reps,
        js_divergences = .result$js_divergences,
        weights = .result$weights,
        global_js_divergence = .result$global_js_divergence,
        p_values = .result$p_values,
        p_values_observed_consensus = .result$p_values_observed_consensus
    )
}
