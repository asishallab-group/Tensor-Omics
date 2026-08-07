# Generated. Do not edit.

#' Grow and track a single ensemble from one seed until a Stop Condition is reached
#'
#' The iteration wrapper described in `misc/mod_STC.md`, "Ensemble identification": each
#' growth step is `grow_ensemble` + `observable`, compared against the last *accepted*
#' iteration via `accept_ensemble` -- skipped, by convention, for the very first growth
#' step, see "First growth step" in the spec -- until one of the four documented Stop
#' Conditions applies. Deliberately sequential: growth at each step depends on the
#' previous one, so there is nothing to parallelize *within* a single seed's growth --
#' outer-level parallelism belongs in whatever calls this once per seed, matching
#' `grow_ensemble_kernel`'s own precedent.
#'
#' Two deliberate readings of the spec, flagged here since the prose leaves them
#' implicit: (1) an isolated seed -- no neighbor at all within its own growth radius --
#' is reported as Stop Condition 4 (a natural fixed point) with a trivial one-member
#' `final_ensemble_mask`, not as "no ensemble" (that is Stop Condition 1's own, distinct
#' meaning). (2) On a rejection (Stop Condition 2 or 3), the *rejected* candidate's
#' observable is still pushed into the trailing history (marked `FALSE` in
#' `accepted_history`) so a caller can see what got rejected and why, even though
#' `final_ensemble_mask` reflects the last *accepted* state, not this one -- otherwise
#' `accepted_history` could only ever read `TRUE`, since only accepted iterations would
#' ever reach the array at all.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering::ensemble_identification}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param seed_index a integer scalar. Index into `vectors`/`kd_indices` of the seed to grow an ensemble around
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param k_min a integer scalar. Neighborhood size for this seed's growth radius, see `calc_ensemble_growth_radius`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#'   The default value is `30`.
#' @param alpha_max a numeric scalar. Maximum tolerated principal angle (radians), see `accept_ensemble`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `2.0 * atan(1.0)`.
#' @param d_max a integer scalar. Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
#'   The minimum valid value is `0`.
#' @param G_max a numeric scalar. Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
#'   The minimum valid value is `0.0`.
#' @param f_max a numeric scalar. Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
#'   The minimum valid value is `above(0.0)`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param a a integer scalar. Minimum accepted-iteration count for a later rejection to count as "stable", see
#'   Stop Condition 2
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param o a integer scalar. Trailing observable-history window depth (`misc/mod_STC.md` suggests 10 as a
#'   sensible default). Always required, never optional with an auto-applied
#'   default here: a Fortran array bound cannot depend on a possibly-absent
#'   optional dummy, and this argument sizes every history output below.
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{final_ensemble_mask}{a logical vector. The last accepted ensemble's membership. All `FALSE` when `stop_reason` is
#'     `STOP_REASON_MAX_SIZE` -- see Stop Condition 1.}
#'   \item{stop_reason}{a integer scalar. Which Stop Condition ended growth: one of
#'     \code{STOP_REASON_MAX_SIZE},
#'     \code{STOP_REASON_REJECTED_AFTER_STABLE},
#'     \code{STOP_REASON_REJECTED_IMMEDIATELY}, or
#'     \code{STOP_REASON_FIXED_POINT} --
#'     or \code{STOP_REASON_ERROR} if
#'     `ierr` is non-zero, in which case every other output for this seed is undefined.}
#'   \item{growth_radius}{a numeric scalar. This seed's growth radius, see `calc_ensemble_growth_radius`}
#'   \item{U_history}{a numeric array of rank 3. Trailing tangent+normal bases, one per retained iteration, oldest to newest;
#'     zero beyond the number of iterations actually retained, see `k_history`}
#'   \item{S_history}{a numeric matrix. Trailing singular values -- not eigenvalues, see "Output" in `misc/mod_STC.md`
#'     -- zero-padded beyond rank and beyond the number of retained iterations}
#'   \item{d_history}{a integer vector. Trailing intrinsic dimensions, one per retained iteration}
#'   \item{G_history}{a numeric vector. Trailing spectral gaps, one per retained iteration}
#'   \item{mu_history}{a numeric matrix. Trailing ensemble centers, one per retained iteration}
#'   \item{k_history}{a integer vector. Trailing ensemble sizes, one per retained iteration. 0 marks a column beyond
#'     the number of iterations actually retained -- a real ensemble size is always
#'     at least 1.}
#'   \item{accepted_history}{a logical vector. Whether the growth iteration retained in the corresponding column was
#'     accepted. Iteration 1 (the bootstrap step) is always `TRUE` by convention.
#'     The single most recent column is `FALSE` when, and only when, growth
#'     stopped via `STOP_REASON_REJECTED_AFTER_STABLE` or
#'     `STOP_REASON_REJECTED_IMMEDIATELY` -- see the module-level note above.}
#'   \item{member_added_at_step}{a integer vector. `MEMBER_ADDED_AT_STEP_NON_MEMBER` for non-members, `MEMBER_ADDED_AT_STEP_SEED`
#'     for the seed itself, the growth-iteration index at which each other member
#'     joined otherwise}
#' @export
ensemble_identification <- function(vectors, kd_indices, dimension_order, seed_index, k_min = 30L, alpha_max, d_max, G_max, f_max = 0.95, a = 2L, o) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    seed_index <- .tox_as_integer_scalar(seed_index, "seed_index")
    k_min <- .tox_as_integer_scalar(k_min, "k_min")
    alpha_max <- .tox_as_double_scalar(alpha_max, "alpha_max")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    f_max <- .tox_as_double_scalar(f_max, "f_max")
    a <- .tox_as_integer_scalar(a, "a")
    o <- .tox_as_integer_scalar(o, "o")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("ensemble_identification_call", vectors, kd_indices, dimension_order, seed_index, k_min, alpha_max, d_max, G_max, f_max, a, o)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_index", "k_min", "alpha_max", "d_max", "G_max", "f_max", "a", "o", "final_ensemble_mask", "stop_reason", "growth_radius", "U_history", "S_history", "d_history", "G_history", "mu_history", "k_history", "accepted_history", "member_added_at_step", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        final_ensemble_mask = .result$final_ensemble_mask,
        stop_reason = .result$stop_reason,
        growth_radius = .result$growth_radius,
        U_history = .result$U_history,
        S_history = .result$S_history,
        d_history = .result$d_history,
        G_history = .result$G_history,
        mu_history = .result$mu_history,
        k_history = .result$k_history,
        accepted_history = .result$accepted_history,
        member_added_at_step = .result$member_added_at_step
    )
}

#' Run ensemble_identification once per seed and assemble the merged, per-ensemble output arrays
#'
#' See `misc/mod_STC.md`, "Ensemble identification", "#### Merged output": every per-seed
#' output of `ensemble_identification` gains one extra trailing dimension of size
#' `n_selected_seed`, one column per seed, in the order seeds occur in `seed_selection_mask`. Each
#' seed's growth is fully independent of every other's, so the per-seed calls below write
#' to disjoint array sections -- safe for `do concurrent`, matching the spec's own "In
#' parallel grow ensembles around each seed vector" and this codebase's existing
#' `do concurrent`-everywhere convention (`tox_shift_vectors_kernel`, `tox_gene_centroids_kernel`,
#' ...). Left for a later pass if it turns out to matter in practice: the spec's own caveat
#' that `do concurrent` is unsafe together with external-library calls under gfortran --
#' `ensemble_identification_kernel` calls LAPACK (`dgesdd`/`dgesvd`) by way of `observable`
#' and `accept_ensemble` -- has not been stress-tested here; `!$omp parallel do` is the
#' documented fallback if it ever is.
#'
#' `ensemble_member_added_at_step` is always collected (see `misc/mod_STC.md`'s "optional,
#' user flag decides" note): unlike Ensemble Reconciliation's JSI, which the user
#' explicitly required to be gated because it adds a real, if small, extra cost per pair,
#' this is just bookkeeping already computed as a side effect of the per-seed growth loop
#' itself -- there is no separate cost left to gate.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering::ensemble_identification_merged}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param seed_selection_mask a logical vector. Seed selection, see `seeds`
#' @param k_min a integer scalar. Neighborhood size for each seed's growth radius, see `calc_ensemble_growth_radius`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors - 1`.
#'   The default value is `30`.
#' @param alpha_max a numeric scalar. Maximum tolerated principal angle (radians), see `accept_ensemble`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `2.0 * atan(1.0)`.
#' @param d_max a integer scalar. Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
#'   The minimum valid value is `0`.
#' @param G_max a numeric scalar. Maximum tolerated |log(G_tp1/G_t)|, see `accept_ensemble`
#'   The minimum valid value is `0.0`.
#' @param f_max a numeric scalar. Ensemble size fraction of N above which growth is abandoned, see Stop Condition 1
#'   The minimum valid value is `above(0.0)`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param a a integer scalar. Minimum accepted-iteration count for a later rejection to count as "stable", see
#'   Stop Condition 2
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param o a integer scalar. Trailing observable-history window depth, see `ensemble_identification`. Always
#'   required, for the same reason as there: it sizes every history output below.
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{ensemble_masks}{a logical matrix. Per-ensemble accepted membership, one column per seed, see `final_ensemble_mask`}
#'   \item{ensemble_stop_reason}{a integer vector. Per-ensemble Stop Condition, see `ensemble_identification`}
#'   \item{ensemble_growth_radii}{a numeric vector. Per-ensemble growth radius, see "Local Radius Identification"}
#'   \item{ensemble_U_history}{a numeric array of rank 4. Per-ensemble trailing tangent+normal bases, see `U_history`}
#'   \item{ensemble_S_history}{a numeric array of rank 3. Per-ensemble trailing singular values, see `S_history`}
#'   \item{ensemble_d_history}{a integer matrix. Per-ensemble trailing intrinsic dimensions, see `d_history`}
#'   \item{ensemble_G_history}{a numeric matrix. Per-ensemble trailing spectral gaps, see `G_history`}
#'   \item{ensemble_mu_history}{a numeric array of rank 3. Per-ensemble trailing centers, see `mu_history`}
#'   \item{ensemble_k_history}{a integer matrix. Per-ensemble trailing sizes, see `k_history`}
#'   \item{ensemble_accepted_history}{a logical matrix. Per-ensemble trailing accepted flags, see `accepted_history`}
#'   \item{ensemble_member_added_at_step}{a integer matrix. Per-ensemble growth-iteration-joined bookkeeping, see `member_added_at_step`}
#' @export
ensemble_identification_merged <- function(vectors, kd_indices, dimension_order, seed_selection_mask, k_min = 30L, alpha_max, d_max, G_max, f_max = 0.95, a = 2L, o) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    seed_selection_mask <- .tox_as_logical(seed_selection_mask, "seed_selection_mask")
    k_min <- .tox_as_integer_scalar(k_min, "k_min")
    alpha_max <- .tox_as_double_scalar(alpha_max, "alpha_max")
    d_max <- .tox_as_integer_scalar(d_max, "d_max")
    G_max <- .tox_as_double_scalar(G_max, "G_max")
    f_max <- .tox_as_double_scalar(f_max, "f_max")
    a <- .tox_as_integer_scalar(a, "a")
    o <- .tox_as_integer_scalar(o, "o")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])
    if (length(seed_selection_mask) != dim(vectors)[2])
        .tox_shape_error("seed_selection_mask", length(seed_selection_mask), "vectors", dim(vectors)[2])

    .result <- .Call("ensemble_identification_merged_call", vectors, kd_indices, dimension_order, seed_selection_mask, k_min, alpha_max, d_max, G_max, f_max, a, o)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "seed_selection_mask", "n_selected_seed", "k_min", "alpha_max", "d_max", "G_max", "f_max", "a", "o", "ensemble_masks", "ensemble_stop_reason", "ensemble_growth_radii", "ensemble_U_history", "ensemble_S_history", "ensemble_d_history", "ensemble_G_history", "ensemble_mu_history", "ensemble_k_history", "ensemble_accepted_history", "ensemble_member_added_at_step", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, "ensemble_masks", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "ensemble_U_history", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        ensemble_masks = .result$ensemble_masks,
        ensemble_stop_reason = .result$ensemble_stop_reason,
        ensemble_growth_radii = .result$ensemble_growth_radii,
        ensemble_U_history = .result$ensemble_U_history,
        ensemble_S_history = .result$ensemble_S_history,
        ensemble_d_history = .result$ensemble_d_history,
        ensemble_G_history = .result$ensemble_G_history,
        ensemble_mu_history = .result$ensemble_mu_history,
        ensemble_k_history = .result$ensemble_k_history,
        ensemble_accepted_history = .result$ensemble_accepted_history,
        ensemble_member_added_at_step = .result$ensemble_member_added_at_step
    )
}
