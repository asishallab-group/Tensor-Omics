# Generated. Do not edit.

#' Identify and group intersecting ensembles from Ensemble Identification's merged output
#'
#' Detecting an intersection at all already requires $|\mathcal{E}_i \cap \mathcal{E}_j|$,
#' needed by every mode; the JSI itself is a single extra $O(1)$ step per pair once each
#' ensemble's own size is known, via
#' $|\mathcal{E}_i \cup \mathcal{E}_j| = |\mathcal{E}_i| + |\mathcal{E}_j| - |\mathcal{E}_i \cap \mathcal{E}_j|$
#' -- but modes 1 and 3 do not need it for their own decision, so its computation is
#' guarded behind `report_jsi`, never unconditional (see `misc/mod_STC.md`'s explicit note
#' on this).
#'
#' Modes 2 and 3 group via a union-find over the qualifying-edge graph (`stc_uf_find`/
#' `stc_uf_union` below), unioning the smaller index under the larger's root so that a
#' component's root is always its own smallest member -- which is what makes the single
#' pass `do r = 1, n_ensembles` below both find every component exactly once and emit them
#' in ascending order of each group's smallest member, with no separate bookkeeping. A
#' discovered component larger than `max_group_size` is a genuine runtime condition no
#' static input check could foresee (it depends on the actual intersection pattern), so it
#' is reported via `ierr` rather than silently truncated -- see `codegen_guide.md` section
#' 5.14.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_reconciliation::ensemble_reconciliation}, whose argument names
#' are the ones an error message reports.
#'
#' @param ensemble_masks a logical matrix. Per-ensemble membership, see Ensemble Identification's merged output
#' @param mode a string, one of "report", "merge_jsi", "merge_any". How intersections are processed
#'
#'   The default value is `1`.
#' @param min_jsi a numeric scalar. Minimum Jaccard Similarity Index for an edge to qualify in mode
#'   \code{'merge_jsi'};
#'   ignored in every other mode
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param report_jsi a logical scalar. Whether to compute and return `super_ensembles_JSI` at all -- see the note
#'   above on this being guarded, not unconditional
#'   The default value is `FALSE`.
#' @param max_group_size a integer scalar. Maximum number of ensembles one super-ensemble (one column of `super_ensembles`)
#'   can hold; sizes its row dimension. `misc/mod_STC.md` suggests
#'   $\min(1024, N_{\mathcal{E}})$ as a sensible default -- always required, never
#'   optional with an auto-applied default here, for the same reason as
#'   `ensemble_identification`'s own `o`: a Fortran array bound cannot depend on a
#'   possibly-absent optional dummy, and a runtime-dependent value like
#'   $\min(1024, N_{\mathcal{E}})$ is not the constant expression an auto-applied
#'   default would need to be either.
#'   The minimum valid value is `2`.
#'   The maximum valid value is `n_ensembles`.
#' @return a named list with elements:
#'   \item{super_ensembles}{a integer matrix. One super-ensemble per column: the 1-indexed column indices of `ensemble_masks`
#'     belonging to that group, padded with 0 (invalid, ensembles are 1-indexed) below
#'     the group's actual size, and 0 in every row of an unused trailing column beyond
#'     `n_super_ensembles`. Sized at $N_{\mathcal{E}}(N_{\mathcal{E}}-1)$, twice mode
#'     \code{'report'}'s
#'     own true worst case ($N_{\mathcal{E}}(N_{\mathcal{E}}-1)/2$, every pair
#'     intersects) -- deliberately not divided by 2: the generator translates this
#'     specification expression close to verbatim into the Python/R bindings, where
#'     `/` on two integers is true division, not Fortran's own truncating integer
#'     division, so a literal `/2` here breaks the generated Python binding (a `float`
#'     where `np.empty`'s shape wants an `int`); see `misc/code_gen_footgun.md`. A
#'     safe, if looser, upper bound for modes 2 and 3 too, whose groups can never
#'     outnumber mode 1's own worst case.}
#'   \item{n_super_ensembles}{a integer scalar. Number of leading columns of `super_ensembles`/`super_ensembles_JSI` actually
#'     filled}
#'   \item{super_ensembles_JSI}{a numeric matrix. Column $l$, row $c_i$: the JSI between the ensembles in `super_ensembles(c_i, l)`
#'     and `super_ensembles(c_i + 1, l)`. All zero unless `report_jsi` was requested --
#'     see the note above.}
#' @export
ensemble_reconciliation <- function(ensemble_masks, mode = "report", min_jsi = 0.1, report_jsi = FALSE, max_group_size) {
    ensemble_masks <- .tox_as_logical(ensemble_masks, "ensemble_masks")
    mode <- .tox_as_mode(mode, "mode", c("report", "merge_jsi", "merge_any"))
    min_jsi <- .tox_as_double_scalar(min_jsi, "min_jsi")
    report_jsi <- .tox_as_logical(report_jsi, "report_jsi")
    max_group_size <- .tox_as_integer_scalar(max_group_size, "max_group_size")
    .result <- .Call("ensemble_reconciliation_call", ensemble_masks, mode, min_jsi, report_jsi, max_group_size)
    .arguments <- c("ensemble_masks", "n_vectors", "n_ensembles", "mode", "min_jsi", "report_jsi", "max_group_size", "super_ensembles", "n_super_ensembles", "super_ensembles_JSI", "ierr")
    .sources <- c(NA_character_, "ensemble_masks", "ensemble_masks", NA_character_, NA_character_, NA_character_, "super_ensembles", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        super_ensembles = .result$super_ensembles,
        n_super_ensembles = .result$n_super_ensembles,
        super_ensembles_JSI = .result$super_ensembles_JSI
    )
}
