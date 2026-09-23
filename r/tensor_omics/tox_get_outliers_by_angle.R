# Generated. Do not edit.

#' Compute each family's mean direction and angular dispersion from expression vectors
#'
#' Every gene's expression vector is scaled to unit length, and a family's mean direction is
#' the normalised sum of its members' unit vectors. With \eqn{R} the length of that sum divided
#' by the number of members, the family's angular dispersion is \eqn{\sigma = \sqrt{-2 \ln R}}:
#' zero when all members point the same way, growing as they spread.
#'
#' The dispersion is computed from the members' own angular deviations -- the angles
#' `compute_angular_deviations` reports -- as \eqn{1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)},
#' which equals the definition because the direction is the members' mean direction, and which
#' loses nothing to cancellation when the members are close.
#'
#' A gene whose expression vector is zero has no direction. It is left out of its family's
#' statistics entirely -- neither counted nor summed. A vector whose components are all
#' subnormal (below `tiny(1.0)`, about `2.2e-308`, in magnitude) counts as zero too, in every
#' build: whether such numbers survive at all depends on the build's floating-point mode. A
#' family whose statistics do not exist gets the zero vector as its direction; see `status`.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_family_direction}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Number of gene families
#' @param expression_vectors a numeric matrix. Expression vector of each gene, one per column
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param min_angular_dispersion a numeric scalar. Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
#'   no dispersion is accepted and every family with statistics is reported.
#'   The default value is `0.0`.
#'   The minimum valid value is `0.0`.
#' @param max_angular_dispersion a numeric scalar. Largest accepted angular dispersion; the default is the dispersion of \eqn{R = 1/2}
#'   The default value is `sqrt(-2.0*log(0.5))`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `5.0`.
#'   At most 5 (\eqn{R \ge e^{-12.5} \approx 4 \times 10^{-6}}): above that, a resultant that cancels only up to rounding would pass as a stable direction.
#' @return a named list with elements:
#'   \item{family_directions}{a numeric matrix. Mean direction of each family as a unit vector, one per column; the zero vector
#'     where the family has fewer than three members with a direction or no stable
#'     direction (see `status`)}
#'   \item{angular_dispersions}{a numeric vector. Angular dispersion \eqn{\sigma = \sqrt{-2 \ln R}} of each family, or
#'     `-1.0` where `status` reports why there is none}
#'   \item{member_counts}{a integer vector. Number of genes of each family that have a direction, i.e. whose expression vector
#'     is not zero -- the members the statistics are taken over}
#'   \item{status}{a integer vector. Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. \code{STAT_TOO_FEW_MEMBERS} where fewer than three genes of the family have a direction: no direction, no dispersion. \code{STAT_NO_STABLE_DIRECTION} where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. \code{STAT_NO_ANGULAR_VARIATION} where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.}
#' @export
compute_family_direction <- function(n_families, expression_vectors, gene_to_fam, min_angular_dispersion = 0.0, max_angular_dispersion = 1.1774100225154747) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    min_angular_dispersion <- .tox_as_double_scalar(min_angular_dispersion, "min_angular_dispersion")
    max_angular_dispersion <- .tox_as_double_scalar(max_angular_dispersion, "max_angular_dispersion")
    if (length(gene_to_fam) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("compute_family_direction_call", n_families, expression_vectors, gene_to_fam, min_angular_dispersion, max_angular_dispersion)
    .arguments <- c("n_axes", "n_genes", "n_families", "expression_vectors", "gene_to_fam", "family_directions", "angular_dispersions", "member_counts", "status", "min_angular_dispersion", "max_angular_dispersion", "ierr")
    .sources <- c("expression_vectors", "expression_vectors", "family_directions", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        family_directions = .result$family_directions,
        angular_dispersions = .result$angular_dispersions,
        member_counts = .result$member_counts,
        status = .result$status
    )
}

#' Compute the angle between each gene's expression vector and its family's mean direction
#'
#' The angle lies in \eqn{[0, \pi]} and does not depend on the length of either vector.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_angular_deviations}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. Expression vector of each gene, one per column
#' @param family_directions a numeric matrix. Mean direction of each family, one per column, as
#'   \code{\link{compute_family_direction}}
#'   gives it; the zero vector marks a family without one
#'   The minimum valid value is `-1.0`.
#'   The maximum valid value is `1.0`.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @return a numeric vector. Angle in radians between each gene and its family's direction, in \eqn{[0, \pi]}, or
#'   `-1.0` where the gene has no family, its expression vector
#'   is zero or all subnormal, or its family has no direction
#' @export
compute_angular_deviations <- function(expression_vectors, family_directions, gene_to_fam) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_directions <- .tox_as_double_matrix(family_directions, "family_directions")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    if (dim(family_directions)[1] != dim(expression_vectors)[1])
        .tox_shape_error("family_directions", dim(family_directions)[1], "expression_vectors", dim(expression_vectors)[1])
    if (length(gene_to_fam) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("compute_angular_deviations_call", expression_vectors, family_directions, gene_to_fam)
    .arguments <- c("n_axes", "n_genes", "n_families", "expression_vectors", "family_directions", "gene_to_fam", "angular_deviations", "ierr")
    .sources <- c("expression_vectors", "expression_vectors", "family_directions", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$angular_deviations
}

#' Compute each family's circular mean angle and angular dispersion from signed angles
#'
#' Every angle is taken as a unit vector \eqn{(\cos\theta, \sin\theta)}. A family's mean angle is
#' the angle of their sum, and with \eqn{R} the length of that sum divided by the number of
#' members, its angular dispersion is \eqn{\sigma = \sqrt{-2 \ln R}}: zero when all members have
#' the same angle, growing as they spread. It is computed from the members' own angular
#' deviations -- the distances `compute_angular_deviations_rap` reports -- as
#' \eqn{1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)}, which loses nothing to cancellation when
#' the members are close.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_family_direction_rap}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Number of gene families
#' @param signed_angles a numeric vector. Signed angle of each gene in radians, in \eqn{(-\pi, \pi]}, as
#'   \code{\link{clock_hand_angles_for_shift_vectors}}
#'   measures it. A gene without an angle is one with no family in `gene_to_fam`; its
#'   entry here is ignored, but must still lie in the range.
#'   The minimum valid value is `above(-PI)`.
#'   The maximum valid value is `PI`.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param min_angular_dispersion a numeric scalar. Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
#'   no dispersion is accepted and every family with statistics is reported.
#'   The default value is `0.0`.
#'   The minimum valid value is `0.0`.
#' @param max_angular_dispersion a numeric scalar. Largest accepted angular dispersion; the default is the dispersion of \eqn{R = 1/2}
#'   The default value is `sqrt(-2.0*log(0.5))`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `5.0`.
#'   At most 5 (\eqn{R \ge e^{-12.5} \approx 4 \times 10^{-6}}): above that, a resultant that cancels only up to rounding would pass as a stable direction.
#' @return a named list with elements:
#'   \item{family_mean_angles}{a numeric vector. Circular mean angle of each family in radians, in \eqn{(-\pi, \pi]}, or
#'     `-4.0` where the family has fewer than three members or no
#'     stable direction (see `status`)}
#'   \item{angular_dispersions}{a numeric vector. Angular dispersion \eqn{\sigma = \sqrt{-2 \ln R}} of each family, or
#'     `-1.0` where `status` reports why there is none}
#'   \item{member_counts}{a integer vector. Number of genes of each family -- the members the statistics are taken over}
#'   \item{status}{a integer vector. Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. \code{STAT_TOO_FEW_MEMBERS} where fewer than three genes belong to the family: no mean angle, no dispersion. \code{STAT_NO_STABLE_DIRECTION} where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. \code{STAT_NO_ANGULAR_VARIATION} where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.}
#' @export
compute_family_direction_rap <- function(n_families, signed_angles, gene_to_fam, min_angular_dispersion = 0.0, max_angular_dispersion = 1.1774100225154747) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    signed_angles <- .tox_as_double_vector(signed_angles, "signed_angles")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    min_angular_dispersion <- .tox_as_double_scalar(min_angular_dispersion, "min_angular_dispersion")
    max_angular_dispersion <- .tox_as_double_scalar(max_angular_dispersion, "max_angular_dispersion")
    if (length(gene_to_fam) != length(signed_angles))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "signed_angles", length(signed_angles))

    .result <- .Call("compute_family_direction_rap_call", n_families, signed_angles, gene_to_fam, min_angular_dispersion, max_angular_dispersion)
    .arguments <- c("n_genes", "n_families", "signed_angles", "gene_to_fam", "family_mean_angles", "angular_dispersions", "member_counts", "status", "min_angular_dispersion", "max_angular_dispersion", "ierr")
    .sources <- c("signed_angles", "family_mean_angles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        family_mean_angles = .result$family_mean_angles,
        angular_dispersions = .result$angular_dispersions,
        member_counts = .result$member_counts,
        status = .result$status
    )
}

#' Compute the absolute angular distance between each gene's signed angle and its family's mean angle
#'
#' The distance is taken around the circle, so it lies in \eqn{[0, \pi]}: an angle just below
#' \eqn{\pi} and one just above \eqn{-\pi} are close. See \code{wrap_angle}.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_angular_deviations_rap}, whose argument names
#' are the ones an error message reports.
#'
#' @param signed_angles a numeric vector. Signed angle of each gene in radians, in \eqn{(-\pi, \pi]}, as
#'   \code{\link{clock_hand_angles_for_shift_vectors}}
#'   measures it. A gene without an angle is one with no family in `gene_to_fam`; its
#'   entry here is ignored, but must still lie in the range.
#'   The minimum valid value is `above(-PI)`.
#'   The maximum valid value is `PI`.
#' @param family_mean_angles a numeric vector. Circular mean angle of each family in radians, as
#'   \code{\link{compute_family_direction_rap}}
#'   gives it, or `-4.0` for a family without one
#'   The minimum valid value is `above(-PI)`.
#'   The maximum valid value is `PI`.
#'   The value `-4.0` is additionally accepted.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @return a numeric vector. Absolute angular distance in radians between each gene's angle and its family's mean
#'   angle, in \eqn{[0, \pi]}, or `-1.0` where the gene has no
#'   family or its family has no mean angle
#' @export
compute_angular_deviations_rap <- function(signed_angles, family_mean_angles, gene_to_fam) {
    signed_angles <- .tox_as_double_vector(signed_angles, "signed_angles")
    family_mean_angles <- .tox_as_double_vector(family_mean_angles, "family_mean_angles")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    if (length(gene_to_fam) != length(signed_angles))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "signed_angles", length(signed_angles))

    .result <- .Call("compute_angular_deviations_rap_call", signed_angles, family_mean_angles, gene_to_fam)
    .arguments <- c("n_genes", "n_families", "signed_angles", "family_mean_angles", "gene_to_fam", "angular_deviations", "ierr")
    .sources <- c("signed_angles", "family_mean_angles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$angular_deviations
}

#' Divide each gene's angular deviation by its family's angular dispersion
#'
#' The relative angular deviation says how far off a gene is in units of its own family's
#' spread, so genes of tight and loose families can be ranked together. It is a plain quotient,
#' with nothing subtracted. It serves both variants: the angular deviations may come from
#' \code{\link{compute_angular_deviations}} or
#' from \code{\link{compute_angular_deviations_rap}}.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_relative_angular_deviations}, whose argument names
#' are the ones an error message reports.
#'
#' @param angular_deviations a numeric vector. Angular deviation of each gene from its family in radians, or
#'   `-1.0` where there is none
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `PI`.
#'   The value `-1.0` is additionally accepted.
#' @param angular_dispersions a numeric vector. Angular dispersion of each family, or `-1.0` where there
#'   is none
#'   The minimum valid value is `0.0`.
#'   The value `-1.0` is additionally accepted.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `angular_deviations`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @return a numeric vector. Angular deviation divided by the family's angular dispersion, at least zero, or
#'   `-1.0` where the gene has no family or no
#'   angular deviation, or its family's dispersion is the sentinel, zero, or below the
#'   smallest normal number `tiny(1.0)` (about `2.2e-308`) -- a dispersion so small that
#'   the quotient could exceed the largest representable number
#' @export
compute_relative_angular_deviations <- function(angular_deviations, angular_dispersions, gene_to_fam) {
    angular_deviations <- .tox_as_double_vector(angular_deviations, "angular_deviations")
    angular_dispersions <- .tox_as_double_vector(angular_dispersions, "angular_dispersions")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    if (length(gene_to_fam) != length(angular_deviations))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "angular_deviations", length(angular_deviations))

    .result <- .Call("compute_relative_angular_deviations_call", angular_deviations, angular_dispersions, gene_to_fam)
    .arguments <- c("n_genes", "n_families", "angular_deviations", "angular_dispersions", "gene_to_fam", "relative_angular_deviations", "ierr")
    .sources <- c("angular_deviations", "angular_dispersions", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$relative_angular_deviations
}

#' Compute the threshold at which a relative angular deviation counts as an outlier
#'
#' The threshold is the `quantile_level` quantile of the relative angular deviations that
#' exist, interpolated linearly between neighbouring values (the "type 7" quantile). The
#' sentinels do not take part.
#'
#' This is a ranking cut, not a statistical test: whatever the data, it picks out the largest
#' values. \code{\link{flag_angle_outliers}}
#' flags only values strictly above the threshold, so with the default level 0.95 at most about
#' the top 5% are flagged. Where many values tie at the quantile position -- a family of
#' identical genes plus one outlier -- the threshold equals the tied value and only the values
#' above it are flagged; where all values are equal, nothing is. With fewer than two values
#' there is nothing to rank, and the threshold is the largest representable number, so nothing
#' is flagged.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::compute_angle_outlier_threshold}, whose argument names
#' are the ones an error message reports.
#'
#' @param relative_angular_deviations a numeric vector. Relative angular deviation of each gene, or `-1.0`
#'   where there is none
#'   The minimum valid value is `0.0`.
#'   The value `-1.0` is additionally accepted.
#' @param quantile_level a numeric scalar. Quantile level in \eqn{[0, 1]} the threshold is taken at
#'   The default value is `0.95`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a numeric scalar. The relative angular deviation a gene must exceed to be an outlier
#' @export
compute_angle_outlier_threshold <- function(relative_angular_deviations, quantile_level = 0.95) {
    relative_angular_deviations <- .tox_as_double_vector(relative_angular_deviations, "relative_angular_deviations")
    quantile_level <- .tox_as_double_scalar(quantile_level, "quantile_level")
    .result <- .Call("compute_angle_outlier_threshold_call", relative_angular_deviations, quantile_level)
    .arguments <- c("n_genes", "relative_angular_deviations", "threshold", "quantile_level", "ierr")
    .sources <- c("relative_angular_deviations", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$threshold
}

#' Flag the genes whose relative angular deviation lies above the threshold
#'
#' A gene is an outlier when its relative angular deviation is positive and strictly greater
#' than `threshold`; a value equal to the threshold is not flagged. So a gene without a value,
#' or with a deviation of zero, is never flagged, whatever the threshold.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::flag_angle_outliers}, whose argument names
#' are the ones an error message reports.
#'
#' @param relative_angular_deviations a numeric vector. Relative angular deviation of each gene, or `-1.0`
#'   where there is none
#'   The minimum valid value is `0.0`.
#'   The value `-1.0` is additionally accepted.
#' @param threshold a numeric scalar. The relative angular deviation a gene must exceed to be an outlier, as
#'   \code{\link{compute_angle_outlier_threshold}}
#'   computes it
#' @return a logical vector. `TRUE` for each gene that is an outlier
#' @export
flag_angle_outliers <- function(relative_angular_deviations, threshold) {
    relative_angular_deviations <- .tox_as_double_vector(relative_angular_deviations, "relative_angular_deviations")
    threshold <- .tox_as_double_scalar(threshold, "threshold")
    .result <- .Call("flag_angle_outliers_call", relative_angular_deviations, threshold)
    .arguments <- c("n_genes", "relative_angular_deviations", "threshold", "is_outlier", "ierr")
    .sources <- c("relative_angular_deviations", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_outlier
}

#' Detect angle outliers among genes from their expression vectors
#'
#' Runs the spherical variant end to end:
#' \code{\link{compute_family_direction}},
#' \code{\link{compute_angular_deviations}},
#' \code{\link{compute_relative_angular_deviations}},
#' \code{\link{compute_angle_outlier_threshold}} and
#' \code{\link{flag_angle_outliers}}.
#'
#' The threshold is a ranking cut, not a statistical test, and only values strictly above it are
#' flagged: where genes tie at the quantile position (identical genes plus one outlier), only
#' the values above the tie are flagged, and where all are equal, none is.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::detect_angle_outliers}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Number of gene families
#' @param expression_vectors a numeric matrix. Expression vector of each gene, one per column
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param quantile_level a numeric scalar. Quantile level in \eqn{[0, 1]} the threshold is taken at
#'   The default value is `0.95`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param min_angular_dispersion a numeric scalar. Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
#'   no dispersion is accepted and every family with statistics is reported.
#'   The default value is `0.0`.
#'   The minimum valid value is `0.0`.
#' @param max_angular_dispersion a numeric scalar. Largest accepted angular dispersion; the default is the dispersion of \eqn{R = 1/2}
#'   The default value is `sqrt(-2.0*log(0.5))`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `5.0`.
#'   At most 5 (\eqn{R \ge e^{-12.5} \approx 4 \times 10^{-6}}): above that, a resultant that cancels only up to rounding would pass as a stable direction.
#' @return a named list with elements:
#'   \item{family_directions}{a numeric matrix. Mean direction of each family as a unit vector, one per column; the zero vector
#'     where the family has fewer than three members with a direction or no stable
#'     direction (see `status`)}
#'   \item{angular_dispersions}{a numeric vector. Angular dispersion \eqn{\sigma = \sqrt{-2 \ln R}} of each family, or
#'     `-1.0` where `status` reports why there is none}
#'   \item{member_counts}{a integer vector. Number of genes of each family that have a direction, i.e. whose expression vector
#'     is not zero -- the members the statistics are taken over}
#'   \item{status}{a integer vector. Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. \code{STAT_TOO_FEW_MEMBERS} where fewer than three genes of the family have a direction: no direction, no dispersion. \code{STAT_NO_STABLE_DIRECTION} where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. \code{STAT_NO_ANGULAR_VARIATION} where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.}
#'   \item{relative_angular_deviations}{a numeric vector. Angle between each gene and its family's direction, divided by the family's angular
#'     dispersion; `-1.0` where `gene_status` is not zero}
#'   \item{threshold}{a numeric scalar. The relative angular deviation a gene must exceed to be an outlier; the largest
#'     representable number when fewer than two genes have a relative angular deviation}
#'   \item{is_outlier}{a logical vector. `TRUE` for each gene whose relative angular deviation is positive and strictly above
#'     `threshold`}
#'   \item{gene_status}{a integer vector. Zero where the gene has a relative angular deviation, otherwise the reason it has
#'     not: \code{STAT_NO_FAMILY} for a gene without a family,
#'     \code{STAT_ZERO_VECTOR} for one whose expression vector is
#'     zero or all subnormal, and otherwise its family's `status`}
#' @export
detect_angle_outliers <- function(n_families, expression_vectors, gene_to_fam, quantile_level = 0.95, min_angular_dispersion = 0.0, max_angular_dispersion = 1.1774100225154747) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    quantile_level <- .tox_as_double_scalar(quantile_level, "quantile_level")
    min_angular_dispersion <- .tox_as_double_scalar(min_angular_dispersion, "min_angular_dispersion")
    max_angular_dispersion <- .tox_as_double_scalar(max_angular_dispersion, "max_angular_dispersion")
    if (length(gene_to_fam) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("detect_angle_outliers_call", n_families, expression_vectors, gene_to_fam, quantile_level, min_angular_dispersion, max_angular_dispersion)
    .arguments <- c("n_axes", "n_genes", "n_families", "expression_vectors", "gene_to_fam", "family_directions", "angular_dispersions", "member_counts", "status", "relative_angular_deviations", "threshold", "is_outlier", "gene_status", "quantile_level", "min_angular_dispersion", "max_angular_dispersion", "ierr")
    .sources <- c("expression_vectors", "expression_vectors", "family_directions", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        family_directions = .result$family_directions,
        angular_dispersions = .result$angular_dispersions,
        member_counts = .result$member_counts,
        status = .result$status,
        relative_angular_deviations = .result$relative_angular_deviations,
        threshold = .result$threshold,
        is_outlier = .result$is_outlier,
        gene_status = .result$gene_status
    )
}

#' Detect angle outliers among genes from their signed angles in a relative axis plane
#'
#' Runs the RAP variant end to end:
#' \code{\link{compute_family_direction_rap}},
#' \code{\link{compute_angular_deviations_rap}},
#' \code{\link{compute_relative_angular_deviations}},
#' \code{\link{compute_angle_outlier_threshold}} and
#' \code{\link{flag_angle_outliers}}.
#'
#' The threshold is a ranking cut, not a statistical test, and only values strictly above it are
#' flagged: where genes tie at the quantile position (identical genes plus one outlier), only
#' the values above the tie are flagged, and where all are equal, none is.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers_by_angle::detect_angle_outliers_rap}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Number of gene families
#' @param signed_angles a numeric vector. Signed angle of each gene in radians, in \eqn{(-\pi, \pi]}, as
#'   \code{\link{clock_hand_angles_for_shift_vectors}}
#'   measures it. A gene without an angle is one with no family in `gene_to_fam`; its
#'   entry here is ignored, but must still lie in the range.
#'   The minimum valid value is `above(-PI)`.
#'   The maximum valid value is `PI`.
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param quantile_level a numeric scalar. Quantile level in \eqn{[0, 1]} the threshold is taken at
#'   The default value is `0.95`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param min_angular_dispersion a numeric scalar. Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
#'   no dispersion is accepted and every family with statistics is reported.
#'   The default value is `0.0`.
#'   The minimum valid value is `0.0`.
#' @param max_angular_dispersion a numeric scalar. Largest accepted angular dispersion; the default is the dispersion of \eqn{R = 1/2}
#'   The default value is `sqrt(-2.0*log(0.5))`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `5.0`.
#'   At most 5 (\eqn{R \ge e^{-12.5} \approx 4 \times 10^{-6}}): above that, a resultant that cancels only up to rounding would pass as a stable direction.
#' @return a named list with elements:
#'   \item{family_mean_angles}{a numeric vector. Circular mean angle of each family in radians, in \eqn{(-\pi, \pi]}, or
#'     `-4.0` where the family has fewer than three members or no
#'     stable direction (see `status`)}
#'   \item{angular_dispersions}{a numeric vector. Angular dispersion \eqn{\sigma = \sqrt{-2 \ln R}} of each family, or
#'     `-1.0` where `status` reports why there is none}
#'   \item{member_counts}{a integer vector. Number of genes of each family -- the members the statistics are taken over}
#'   \item{status}{a integer vector. Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. \code{STAT_TOO_FEW_MEMBERS} where fewer than three genes belong to the family: no mean angle, no dispersion. \code{STAT_NO_STABLE_DIRECTION} where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. \code{STAT_NO_ANGULAR_VARIATION} where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.}
#'   \item{relative_angular_deviations}{a numeric vector. Absolute angular distance between each gene's angle and its family's mean angle,
#'     divided by the family's angular dispersion; `-1.0`
#'     where `gene_status` is not zero}
#'   \item{threshold}{a numeric scalar. The relative angular deviation a gene must exceed to be an outlier; the largest
#'     representable number when fewer than two genes have a relative angular deviation}
#'   \item{is_outlier}{a logical vector. `TRUE` for each gene whose relative angular deviation is positive and strictly above
#'     `threshold`}
#'   \item{gene_status}{a integer vector. Zero where the gene has a relative angular deviation, otherwise the reason it has
#'     not: \code{STAT_NO_FAMILY} for a gene without a family, and
#'     otherwise its family's `status`}
#' @export
detect_angle_outliers_rap <- function(n_families, signed_angles, gene_to_fam, quantile_level = 0.95, min_angular_dispersion = 0.0, max_angular_dispersion = 1.1774100225154747) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    signed_angles <- .tox_as_double_vector(signed_angles, "signed_angles")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    quantile_level <- .tox_as_double_scalar(quantile_level, "quantile_level")
    min_angular_dispersion <- .tox_as_double_scalar(min_angular_dispersion, "min_angular_dispersion")
    max_angular_dispersion <- .tox_as_double_scalar(max_angular_dispersion, "max_angular_dispersion")
    if (length(gene_to_fam) != length(signed_angles))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "signed_angles", length(signed_angles))

    .result <- .Call("detect_angle_outliers_rap_call", n_families, signed_angles, gene_to_fam, quantile_level, min_angular_dispersion, max_angular_dispersion)
    .arguments <- c("n_genes", "n_families", "signed_angles", "gene_to_fam", "family_mean_angles", "angular_dispersions", "member_counts", "status", "relative_angular_deviations", "threshold", "is_outlier", "gene_status", "quantile_level", "min_angular_dispersion", "max_angular_dispersion", "ierr")
    .sources <- c("signed_angles", "family_mean_angles", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        family_mean_angles = .result$family_mean_angles,
        angular_dispersions = .result$angular_dispersions,
        member_counts = .result$member_counts,
        status = .result$status,
        relative_angular_deviations = .result$relative_angular_deviations,
        threshold = .result$threshold,
        is_outlier = .result$is_outlier,
        gene_status = .result$gene_status
    )
}
