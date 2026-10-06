#include <src/macros.h>

!> Gene outliers, from how far each gene's direction deviates from its family's.
!|
!| Two variants answer the same question on different data. The spherical variant takes
!| expression vectors: every gene points somewhere in expression space, its family has a mean
!| direction, and a gene's angular deviation is the angle between the two. The RAP variant takes
!| one signed angle per gene within a relative axis plane, as
!| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
!| measures them: the family has a circular mean angle, and a gene's angular deviation is its
!| absolute wrapped distance to that mean.
!|
!| Both then share the rest. A family's spread is its angular dispersion
!| \(\sigma = \sqrt{-2 \ln R}\), where \(R\) is the mean resultant length of its members' unit
!| vectors. A gene's relative angular deviation is its angular deviation divided by its family's
!| angular dispersion, so a gene in a tight family counts as far off sooner than one in a loose
!| family. The genes whose relative angular deviation lies strictly above a quantile of all of
!| them are flagged: a ranking cut, not a statistical test, and values equal to the quantile are
!| not flagged.
!|
!| `detect_angle_outliers` and `detect_angle_outliers_rap` run the whole method in one call and
!| are the entry points to reach for first. Each step is callable on its own too:
!| `compute_family_direction` (or `compute_family_direction_rap`), `compute_angular_deviations`
!| (or `compute_angular_deviations_rap`), then the shared `compute_relative_angular_deviations`,
!| `compute_angle_outlier_threshold` and `flag_angle_outliers`.
!|
!| A family's statistics exist only when at least three of its genes have a direction and their
!| unit vectors do not cancel out. A family for which that fails, or whose dispersion lies
!| outside the accepted bounds, is reported in `status` (without a direction unless its only fault
!| is too little angular variation) with a `STAT_*` code from
!| [[tox_errors(module)]], and its genes get no relative angular deviation and are never flagged.
!| A value that does not exist is marked with a sentinel outside the range of real ones: angular
!| deviations, dispersions and relative angular deviations are never negative and use `-1.0`;
!| signed angles and family mean angles lie in \((-\pi, \pi]\), and a missing mean angle is
!| `-4.0`.
!|
!| A dispersion too small to tell from rounding counts as no angular variation, even with a
!| minimum of 0: identical or proportional members compute to deviations of rounding size, not
!| to exactly zero. The limit is about \((m + 1.5\,n + 16) \cdot 2.2 \times 10^{-16}\) rad for
!| a family of \(m\) members in \(n\) axes (no axis term for signed angles) -- about
!| \(10^{-14}\) rad for small families, still only about \(2 \times 10^{-11}\) rad for \(10^5\) members.
module tox_get_outliers_by_angle_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, is_err, STAT_NO_STABLE_DIRECTION, STAT_NO_ANGULAR_VARIATION, &
                          STAT_TOO_FEW_MEMBERS, STAT_NO_FAMILY, STAT_ZERO_VECTOR
    use f42_math_impl, only: wrap_angle, angular_distance, one_minus_cosine, log_one_minus
    use f42_vector_impl, only: scaled_length, has_direction, angle_to_direction
    use f42_sort_impl, only: sort_array_heapsort
    use f42_stats_impl, only: calc_percentile_impl
    M_IMPLICIT_NONE

    private
    public :: compute_family_direction_impl, compute_angular_deviations_impl, &
              compute_family_direction_rap_impl, compute_angular_deviations_rap_impl, &
              compute_relative_angular_deviations_impl, compute_angle_outlier_threshold_impl, &
              flag_angle_outliers_impl, detect_angle_outliers_impl, detect_angle_outliers_rap_impl

! Every sentinel is negative, the impossible range of the value it stands in for. The mean angle
! of a family is an angle in (-pi, pi], so its sentinel lies below -pi instead.
#define CM_ANGULAR_DEVIATIONS_SENTINEL -1.0_real64
#define CM_ANGULAR_DISPERSION_SENTINEL -1.0_real64
#define CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL -1.0_real64
#define CM_FAMILY_MEAN_ANGLE_SENTINEL -4.0_real64

! The accepted angular dispersions: any above the rounding noise floor, at most the dispersion of
! R = 1/2. The minimum defaults to 0 rather than a floor like 0.01, because a fixed floor hides
! real outliers: 999 identical genes plus one gene 0.3 rad off have a dispersion of 9.45e-3, so a
! minimum of 0.01 rejected the family and the one clear outlier was never flagged. Identical
! members are caught by the noise floor instead (see dispersion_noise_floor).
#define CM_MIN_ANGULAR_DISPERSION_DEFAULT 0.0_real64
#define CM_MAX_ANGULAR_DISPERSION_DEFAULT sqrt(-2.0_real64*log(0.5_real64))
! The largest maximum a caller may set, see max_angular_dispersion
#define CM_MAX_ANGULAR_DISPERSION_LIMIT 5.0_real64

#define CM_OUTLIER_QUANTILE_LEVEL_DEFAULT 0.95_real64

! A family needs this many genes with a direction before it has statistics
#define CM_MIN_FAMILY_MEMBERS 3_int32

! The `status` of a family, for both variants: what the family lacks is DIRECTION
#define CM_FAMILY_STATUS_DOC(DIRECTION, TOO_FEW) Zero where the family has a DIRECTION and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where TOO_FEW: no DIRECTION, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no DIRECTION, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the DIRECTION kept.

contains

    !> summary: Compute each family's mean direction and angular dispersion from expression vectors
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Every gene's expression vector is scaled to unit length, and a family's mean direction is
    !| the normalised sum of its members' unit vectors. With \(R\) the length of that sum divided
    !| by the number of members, the family's angular dispersion is \(\sigma = \sqrt{-2 \ln R}\):
    !| zero when all members point the same way, growing as they spread.
    !|
    !| The dispersion is computed from the members' own angular deviations -- the angles
    !| `compute_angular_deviations` reports -- as \(1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)\),
    !| which equals the definition because the direction is the members' mean direction, and which
    !| loses nothing to cancellation when the members are close.
    !|
    !| A gene whose expression vector is zero has no direction. It is left out of its family's
    !| statistics entirely -- neither counted nor summed. A vector whose components are all
    !| subnormal (below `tiny(1.0)`, about `2.2e-308`, in magnitude) counts as zero too, in every
    !| build: whether such numbers survive at all depends on the build's floating-point mode. A
    !| family whose statistics do not exist gets the zero vector as its direction; see `status`.
    pure subroutine compute_family_direction_impl(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                                  family_directions, angular_dispersions, member_counts, status, &
                                                  min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
            !! Expression vector of each gene, one per column
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: family_directions(n_axes, n_families)
            !! Mean direction of each family as a unit vector, one per column; the zero vector
            !! where the family has fewer than three members with a direction or no stable
            !! direction (see `status`)
        real(real64), intent(out) :: angular_dispersions(n_families)
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `CM_ANGULAR_DISPERSION_SENTINEL` where `status` reports why there is none
        integer(int32), intent(out) :: member_counts(n_families)
            !! Number of genes of each family that have a direction, i.e. whose expression vector
            !! is not zero -- the members the statistics are taken over
        integer(int32), intent(out) :: status(n_families)
            !! CM_FAMILY_STATUS_DOC(direction, fewer than three genes of the family have a direction)
        real(real64), intent(in), optional :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! DM_DEFAULT(CM_MIN_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! DM_DEFAULT(CM_MAX_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(CM_MAX_ANGULAR_DISPERSION_LIMIT)
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

        real(real64) :: actual_min_dispersion, actual_max_dispersion
        real(real64) :: gene_factor, gene_scaled_norm, sum_factor, sum_scaled_norm, deviation
        integer(int32) :: i_gene, i_family, i_axis

        M_DEFAULT_VAL(min_angular_dispersion, actual_min_dispersion, CM_MIN_ANGULAR_DISPERSION_DEFAULT)
        M_DEFAULT_VAL(max_angular_dispersion, actual_max_dispersion, CM_MAX_ANGULAR_DISPERSION_DEFAULT)

        family_directions = 0.0_real64
        member_counts = 0_int32

        ! scatters into the families, so it stays sequential over the genes
        do i_gene = 1, n_genes
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle

            call scaled_length(n_axes, expression_vectors(:, i_gene), gene_factor, gene_scaled_norm)
            if (gene_scaled_norm == 0.0_real64) cycle  ! a zero vector has no direction

            member_counts(i_family) = member_counts(i_family) + 1_int32
            do concurrent (i_axis = 1:n_axes) shared(family_directions, expression_vectors, i_family, i_gene, &
                                                      gene_factor, gene_scaled_norm)
                family_directions(i_axis, i_family) = family_directions(i_axis, i_family) &
                                                      + (expression_vectors(i_axis, i_gene)*gene_factor)/gene_scaled_norm
            end do
        end do

        ! the mean direction of every family with enough members. Each iteration reads and writes
        ! only its own family's column and entries; `angular_dispersions` then gathers the members'
        ! 1 - cos(deviation) until the dispersions replace it
        do concurrent (i_family = 1:n_families) local(sum_factor, sum_scaled_norm) &
            shared(family_directions, angular_dispersions, member_counts, status, n_axes)

            angular_dispersions(i_family) = 0.0_real64
            call set_ok(status(i_family))
            if (member_counts(i_family) < CM_MIN_FAMILY_MEMBERS) then
                family_directions(:, i_family) = 0.0_real64
                call set_err(status(i_family), STAT_TOO_FEW_MEMBERS)
                cycle
            end if

            call scaled_length(n_axes, family_directions(:, i_family), sum_factor, sum_scaled_norm)
            if (sum_scaled_norm == 0.0_real64) then  ! the unit vectors cancel out exactly
                family_directions(:, i_family) = 0.0_real64
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if
            family_directions(:, i_family) = (family_directions(:, i_family)*sum_factor)/sum_scaled_norm
        end do

        ! each member's deviation from its family's direction, the very angle
        ! compute_angular_deviations reports; scatters, so it stays sequential
        do i_gene = 1, n_genes
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle
            if (is_err(status(i_family))) cycle

            deviation = angle_to_direction(n_axes, expression_vectors(:, i_gene), family_directions(:, i_family))
            if (deviation < 0.0_real64) cycle  ! a zero vector has no direction
            angular_dispersions(i_family) = angular_dispersions(i_family) + one_minus_cosine(deviation)
        end do

        call classify_dispersions(n_families, member_counts, n_axes, actual_min_dispersion, actual_max_dispersion, &
                                  angular_dispersions, status)
        ! a direction that is not stable is no direction: the mean of a nearly cancelling
        ! resultant is rounding noise
        do concurrent (i_family = 1:n_families) shared(family_directions, status)
            if (status(i_family) == STAT_NO_STABLE_DIRECTION) family_directions(:, i_family) = 0.0_real64
        end do
    end subroutine compute_family_direction_impl

    !> summary: Compute the angle between each gene's expression vector and its family's mean direction
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The angle lies in \([0, \pi]\) and does not depend on the length of either vector.
    pure subroutine compute_angular_deviations_impl(n_axes, n_genes, n_families, expression_vectors, family_directions, &
                                                    gene_to_fam, angular_deviations)
        integer(int32), intent(in) :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
            !! Expression vector of each gene, one per column
        real(real64), intent(in) :: family_directions(n_axes, n_families)
            !! Mean direction of each family, one per column, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_family_direction_impl(subroutine)]]
            !! gives it; the zero vector marks a family without one
            !! DM_MIN(-1.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: angular_deviations(n_genes)
            !! Angle in radians between each gene and its family's direction, in \([0, \pi]\), or
            !! `CM_ANGULAR_DEVIATIONS_SENTINEL` where the gene has no family, its expression vector
            !! is zero or all subnormal, or its family has no direction

        real(real64) :: deviation
        integer(int32) :: i_gene, i_family

        do concurrent (i_gene = 1:n_genes) local(i_family, deviation) &
            shared(expression_vectors, family_directions, gene_to_fam, angular_deviations, n_axes)

            angular_deviations(i_gene) = CM_ANGULAR_DEVIATIONS_SENTINEL
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle

            deviation = angle_to_direction(n_axes, expression_vectors(:, i_gene), family_directions(:, i_family))
            if (deviation < 0.0_real64) cycle  ! either vector is zero: no direction
            angular_deviations(i_gene) = deviation
        end do
    end subroutine compute_angular_deviations_impl

    !> summary: Compute each family's circular mean angle and angular dispersion from signed angles
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Every angle is taken as a unit vector \((\cos\theta, \sin\theta)\). A family's mean angle is
    !| the angle of their sum, and with \(R\) the length of that sum divided by the number of
    !| members, its angular dispersion is \(\sigma = \sqrt{-2 \ln R}\): zero when all members have
    !| the same angle, growing as they spread. It is computed from the members' own angular
    !| deviations -- the distances `compute_angular_deviations_rap` reports -- as
    !| \(1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)\), which loses nothing to cancellation when
    !| the members are close.
    pure subroutine compute_family_direction_rap_impl(n_genes, n_families, signed_angles, gene_to_fam, &
                                                      family_mean_angles, angular_dispersions, member_counts, status, &
                                                      min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: signed_angles(n_genes)
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! DM_MIN(above(-PI))
            !! DM_MAX(PI)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(signed_angles)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: family_mean_angles(n_families)
            !! Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            !! `CM_FAMILY_MEAN_ANGLE_SENTINEL` where the family has fewer than three members or no
            !! stable direction (see `status`)
        real(real64), intent(out) :: angular_dispersions(n_families)
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `CM_ANGULAR_DISPERSION_SENTINEL` where `status` reports why there is none
        integer(int32), intent(out) :: member_counts(n_families)
            !! Number of genes of each family -- the members the statistics are taken over
        integer(int32), intent(out) :: status(n_families)
            !! CM_FAMILY_STATUS_DOC(mean angle, fewer than three genes belong to the family)
        real(real64), intent(in), optional :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! DM_DEFAULT(CM_MIN_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! DM_DEFAULT(CM_MAX_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(CM_MAX_ANGULAR_DISPERSION_LIMIT)
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

        real(real64) :: actual_min_dispersion, actual_max_dispersion
        real(real64) :: sum_sin, sum_cos
        integer(int32) :: i_gene, i_family

        M_DEFAULT_VAL(min_angular_dispersion, actual_min_dispersion, CM_MIN_ANGULAR_DISPERSION_DEFAULT)
        M_DEFAULT_VAL(max_angular_dispersion, actual_max_dispersion, CM_MAX_ANGULAR_DISPERSION_DEFAULT)

        ! the two outputs hold the sums of sines and cosines until the loop over the families
        ! replaces them with the statistics
        family_mean_angles = 0.0_real64
        angular_dispersions = 0.0_real64
        member_counts = 0_int32

        ! scatters into the families, so it stays sequential over the genes
        do i_gene = 1, n_genes
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle

            member_counts(i_family) = member_counts(i_family) + 1_int32
            family_mean_angles(i_family) = family_mean_angles(i_family) + sin(signed_angles(i_gene))
            angular_dispersions(i_family) = angular_dispersions(i_family) + cos(signed_angles(i_gene))
        end do

        ! each iteration reads its own family's sums before it overwrites them; `angular_dispersions`
        ! then gathers the members' 1 - cos(deviation) until the dispersions replace it
        do concurrent (i_family = 1:n_families) local(sum_sin, sum_cos) &
            shared(family_mean_angles, angular_dispersions, member_counts, status)

            sum_sin = family_mean_angles(i_family)
            sum_cos = angular_dispersions(i_family)
            angular_dispersions(i_family) = 0.0_real64
            call set_ok(status(i_family))

            if (member_counts(i_family) < CM_MIN_FAMILY_MEMBERS) then
                family_mean_angles(i_family) = CM_FAMILY_MEAN_ANGLE_SENTINEL
                call set_err(status(i_family), STAT_TOO_FEW_MEMBERS)
                cycle
            end if

            ! Library sines and cosines cancel to exactly zero only for special families (such
            ! as x, -x, pi - x, x - pi for some x), and such a resultant has no angle at all
            if (hypot(sum_sin, sum_cos) == 0.0_real64) then
                family_mean_angles(i_family) = CM_FAMILY_MEAN_ANGLE_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if

            ! atan2 returns exactly -pi where the sine sum is a negative number too small to move the
            ! result off -pi (members at pi and just above -pi); the same direction is pi
            family_mean_angles(i_family) = wrap_angle(atan2(sum_sin, sum_cos))
        end do

        ! each member's deviation from its family's mean, the very distance
        ! compute_angular_deviations_rap reports; scatters, so it stays sequential
        do i_gene = 1, n_genes
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle
            if (is_err(status(i_family))) cycle

            angular_dispersions(i_family) = angular_dispersions(i_family) &
                + one_minus_cosine(angular_distance(signed_angles(i_gene), family_mean_angles(i_family)))
        end do

        ! signed angles have no axes, so the noise floor has no axis term
        call classify_dispersions(n_families, member_counts, 0_int32, actual_min_dispersion, actual_max_dispersion, &
                                  angular_dispersions, status)
        ! a mean that is not stable is no mean: that of a nearly cancelling resultant is noise
        do concurrent (i_family = 1:n_families) shared(family_mean_angles, status)
            if (status(i_family) == STAT_NO_STABLE_DIRECTION) family_mean_angles(i_family) = CM_FAMILY_MEAN_ANGLE_SENTINEL
        end do
    end subroutine compute_family_direction_rap_impl

    !> summary: Compute the absolute angular distance between each gene's signed angle and its family's mean angle
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The distance is taken around the circle, so it lies in \([0, \pi]\): an angle just below
    !| \(\pi\) and one just above \(-\pi\) are close. See [[f42_math_impl(module):wrap_angle(function)]].
    pure subroutine compute_angular_deviations_rap_impl(n_genes, n_families, signed_angles, family_mean_angles, &
                                                        gene_to_fam, angular_deviations)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: signed_angles(n_genes)
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! DM_MIN(above(-PI))
            !! DM_MAX(PI)
        real(real64), intent(in) :: family_mean_angles(n_families)
            !! Circular mean angle of each family in radians, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_family_direction_rap_impl(subroutine)]]
            !! gives it, or `CM_FAMILY_MEAN_ANGLE_SENTINEL` for a family without one
            !! DM_MIN(above(-PI))
            !! DM_MAX(PI)
            !! DM_SENTINEL(CM_FAMILY_MEAN_ANGLE_SENTINEL)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(signed_angles)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: angular_deviations(n_genes)
            !! Absolute angular distance in radians between each gene's angle and its family's mean
            !! angle, in \([0, \pi]\), or `CM_ANGULAR_DEVIATIONS_SENTINEL` where the gene has no
            !! family or its family has no mean angle

        integer(int32) :: i_gene, i_family

        do concurrent (i_gene = 1:n_genes) local(i_family) &
            shared(signed_angles, family_mean_angles, gene_to_fam, angular_deviations)

            angular_deviations(i_gene) = CM_ANGULAR_DEVIATIONS_SENTINEL
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle
            if (family_mean_angles(i_family) == CM_FAMILY_MEAN_ANGLE_SENTINEL) cycle

            angular_deviations(i_gene) = angular_distance(signed_angles(i_gene), family_mean_angles(i_family))
        end do
    end subroutine compute_angular_deviations_rap_impl

    !> summary: Divide each gene's angular deviation by its family's angular dispersion
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The relative angular deviation says how far off a gene is in units of its own family's
    !| spread, so genes of tight and loose families can be ranked together. It is a plain quotient,
    !| with nothing subtracted. It serves both variants: the angular deviations may come from
    !| [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_impl(subroutine)]] or
    !| from [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_rap_impl(subroutine)]].
    pure subroutine compute_relative_angular_deviations_impl(n_genes, n_families, angular_deviations, angular_dispersions, &
                                                             gene_to_fam, relative_angular_deviations)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: angular_deviations(n_genes)
            !! Angular deviation of each gene from its family in radians, or
            !! `CM_ANGULAR_DEVIATIONS_SENTINEL` where there is none
            !! DM_MIN(0.0_real64)
            !! DM_MAX(PI)
            !! DM_SENTINEL(CM_ANGULAR_DEVIATIONS_SENTINEL)
        real(real64), intent(in) :: angular_dispersions(n_families)
            !! Angular dispersion of each family, or `CM_ANGULAR_DISPERSION_SENTINEL` where there
            !! is none
            !! DM_MIN(0.0_real64)
            !! DM_SENTINEL(CM_ANGULAR_DISPERSION_SENTINEL)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(angular_deviations)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: relative_angular_deviations(n_genes)
            !! Angular deviation divided by the family's angular dispersion, at least zero, or
            !! `CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL` where the gene has no family or no
            !! angular deviation, or its family's dispersion is the sentinel, zero, or below the
            !! smallest normal number `tiny(1.0)` (about `2.2e-308`) -- a dispersion so small that
            !! the quotient could exceed the largest representable number

        real(real64) :: dispersion
        integer(int32) :: i_gene, i_family

        do concurrent (i_gene = 1:n_genes) local(i_family, dispersion) &
            shared(angular_deviations, angular_dispersions, gene_to_fam, relative_angular_deviations)

            relative_angular_deviations(i_gene) = CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL
            i_family = gene_to_fam(i_gene)
            if (i_family == M_GENE_TO_FAM_SENTINEL) cycle
            if (angular_deviations(i_gene) < 0.0_real64) cycle

            ! an angle is at most pi, so from `tiny` on the quotient stays below `huge`; the
            ! sentinel and zero fall below it too
            dispersion = angular_dispersions(i_family)
            if (dispersion < tiny(1.0_real64)) cycle

            relative_angular_deviations(i_gene) = angular_deviations(i_gene)/dispersion
        end do
    end subroutine compute_relative_angular_deviations_impl

    !> summary: Compute the threshold at which a relative angular deviation counts as an outlier
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The threshold is the `quantile_level` quantile of the relative angular deviations that
    !| exist, interpolated linearly between neighbouring values (the "type 7" quantile). The
    !| sentinels do not take part.
    !|
    !| This is a ranking cut, not a statistical test: whatever the data, it picks out the largest
    !| values. [[tox_get_outliers_by_angle_impl(module):flag_angle_outliers_impl(subroutine)]]
    !| flags only values strictly above the threshold, so with the default level 0.95 at most about
    !| the top 5% are flagged. Where many values tie at the quantile position -- a family of
    !| identical genes plus one outlier -- the threshold equals the tied value and only the values
    !| above it are flagged; where all values are equal, nothing is. With fewer than two values
    !| there is nothing to rank, and the threshold is the largest representable number, so nothing
    !| is flagged.
    pure subroutine compute_angle_outlier_threshold_impl(n_genes, relative_angular_deviations, &
                                                         tmp_relative_angular_deviations_perm, threshold, quantile_level)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), intent(in) :: relative_angular_deviations(n_genes)
            !! Relative angular deviation of each gene, or `CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL`
            !! where there is none
            !! DM_MIN(0.0_real64)
            !! DM_SENTINEL(CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL)
        integer(int32), intent(out) :: tmp_relative_angular_deviations_perm(n_genes)
            !! Work array: the indices of the values that exist, sorted by value
        real(real64), intent(out) :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier
        real(real64), intent(in), optional :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! DM_DEFAULT(CM_OUTLIER_QUANTILE_LEVEL_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)

        real(real64) :: actual_quantile_level
        integer(int32) :: i_gene, n_valid

        M_DEFAULT_VAL(quantile_level, actual_quantile_level, CM_OUTLIER_QUANTILE_LEVEL_DEFAULT)

        ! gather the indices of the values that exist to the front, so the quantile can be taken
        ! over those leading entries alone
        n_valid = 0_int32
        do i_gene = 1, n_genes
            if (relative_angular_deviations(i_gene) >= 0.0_real64) then
                n_valid = n_valid + 1_int32
                tmp_relative_angular_deviations_perm(n_valid) = i_gene
            end if
        end do
        tmp_relative_angular_deviations_perm(n_valid + 1:) = 0_int32

        if (n_valid < 2_int32) then
            threshold = huge(1.0_real64)
            return
        end if

        call sort_array_heapsort(relative_angular_deviations, tmp_relative_angular_deviations_perm(1:n_valid))
        call calc_percentile_impl(relative_angular_deviations, n_genes, tmp_relative_angular_deviations_perm, &
                                  actual_quantile_level, threshold, n_considered=n_valid)
    end subroutine compute_angle_outlier_threshold_impl

    !> summary: Flag the genes whose relative angular deviation lies above the threshold
    !| AUTHOR_FRANZ_ERIC_SILL
    !| A gene is an outlier when its relative angular deviation is positive and strictly greater
    !| than `threshold`; a value equal to the threshold is not flagged. So a gene without a value,
    !| or with a deviation of zero, is never flagged, whatever the threshold.
    pure subroutine flag_angle_outliers_impl(n_genes, relative_angular_deviations, threshold, is_outlier)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), intent(in) :: relative_angular_deviations(n_genes)
            !! Relative angular deviation of each gene, or `CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL`
            !! where there is none
            !! DM_MIN(0.0_real64)
            !! DM_SENTINEL(CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL)
        real(real64), intent(in) :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_angle_outlier_threshold_impl(subroutine)]]
            !! computes it
        logical(c_bool), intent(out) :: is_outlier(n_genes)
            !! `.true.` for each gene that is an outlier

        integer(int32) :: i_gene

        ! the flag rule lives in this one statement: strictly above, so values tied at the
        ! threshold are not flagged
        do concurrent (i_gene = 1:n_genes) shared(relative_angular_deviations, threshold, is_outlier)
            is_outlier(i_gene) = relative_angular_deviations(i_gene) > 0.0_real64 &
                                 .and. relative_angular_deviations(i_gene) > threshold
        end do
    end subroutine flag_angle_outliers_impl

    !> summary: Detect angle outliers among genes from their expression vectors
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Runs the spherical variant end to end:
    !| [[tox_get_outliers_by_angle_impl(module):compute_family_direction_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_relative_angular_deviations_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_angle_outlier_threshold_impl(subroutine)]] and
    !| [[tox_get_outliers_by_angle_impl(module):flag_angle_outliers_impl(subroutine)]].
    !|
    !| The threshold is a ranking cut, not a statistical test, and only values strictly above it are
    !| flagged: where genes tie at the quantile position (identical genes plus one outlier), only
    !| the values above the tie are flagged, and where all are equal, none is.
    pure subroutine detect_angle_outliers_impl(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                               family_directions, angular_dispersions, member_counts, status, &
                                               relative_angular_deviations, threshold, is_outlier, gene_status, &
                                               tmp_angular_deviations, tmp_relative_angular_deviations_perm, &
                                               quantile_level, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: expression_vectors(n_axes, n_genes)
            !! Expression vector of each gene, one per column
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(expression_vectors)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: family_directions(n_axes, n_families)
            !! Mean direction of each family as a unit vector, one per column; the zero vector
            !! where the family has fewer than three members with a direction or no stable
            !! direction (see `status`)
        real(real64), intent(out) :: angular_dispersions(n_families)
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `CM_ANGULAR_DISPERSION_SENTINEL` where `status` reports why there is none
        integer(int32), intent(out) :: member_counts(n_families)
            !! Number of genes of each family that have a direction, i.e. whose expression vector
            !! is not zero -- the members the statistics are taken over
        integer(int32), intent(out) :: status(n_families)
            !! CM_FAMILY_STATUS_DOC(direction, fewer than three genes of the family have a direction)
        real(real64), intent(out) :: relative_angular_deviations(n_genes)
            !! Angle between each gene and its family's direction, divided by the family's angular
            !! dispersion; `CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL` where `gene_status` is not zero
        real(real64), intent(out) :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), intent(out) :: is_outlier(n_genes)
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(int32), intent(out) :: gene_status(n_genes)
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family,
            !! [[tox_errors(module):STAT_ZERO_VECTOR(variable)]] for one whose expression vector is
            !! zero or all subnormal, and otherwise its family's `status`
        real(real64), intent(out) :: tmp_angular_deviations(n_genes)
            !! Work array: the angle between each gene and its family's direction
        integer(int32), intent(out) :: tmp_relative_angular_deviations_perm(n_genes)
            !! Work array: the indices of the relative angular deviations that exist, sorted by value
        real(real64), intent(in), optional :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! DM_DEFAULT(CM_OUTLIER_QUANTILE_LEVEL_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(in), optional :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! DM_DEFAULT(CM_MIN_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! DM_DEFAULT(CM_MAX_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(CM_MAX_ANGULAR_DISPERSION_LIMIT)
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

        integer(int32) :: i_gene

        call compute_family_direction_impl(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                           family_directions, angular_dispersions, member_counts, status, &
                                           min_angular_dispersion, max_angular_dispersion)
        call compute_angular_deviations_impl(n_axes, n_genes, n_families, expression_vectors, family_directions, &
                                             gene_to_fam, tmp_angular_deviations)
        call compute_relative_angular_deviations_impl(n_genes, n_families, tmp_angular_deviations, angular_dispersions, &
                                                      gene_to_fam, relative_angular_deviations)
        call compute_angle_outlier_threshold_impl(n_genes, relative_angular_deviations, &
                                                  tmp_relative_angular_deviations_perm, threshold, quantile_level)
        call flag_angle_outliers_impl(n_genes, relative_angular_deviations, threshold, is_outlier)

        do concurrent (i_gene = 1:n_genes) shared(gene_status, gene_to_fam, status, expression_vectors)
            gene_status(i_gene) = gene_status_of(gene_to_fam(i_gene), status, has_direction(n_axes, expression_vectors(:, i_gene)))
        end do
    end subroutine detect_angle_outliers_impl

    !> summary: Detect angle outliers among genes from their signed angles in a relative axis plane
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Runs the RAP variant end to end:
    !| [[tox_get_outliers_by_angle_impl(module):compute_family_direction_rap_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_rap_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_relative_angular_deviations_impl(subroutine)]],
    !| [[tox_get_outliers_by_angle_impl(module):compute_angle_outlier_threshold_impl(subroutine)]] and
    !| [[tox_get_outliers_by_angle_impl(module):flag_angle_outliers_impl(subroutine)]].
    !|
    !| The threshold is a ranking cut, not a statistical test, and only values strictly above it are
    !| flagged: where genes tie at the quantile position (identical genes plus one outlier), only
    !| the values above the tie are flagged, and where all are equal, none is.
    pure subroutine detect_angle_outliers_rap_impl(n_genes, n_families, signed_angles, gene_to_fam, &
                                                   family_mean_angles, angular_dispersions, member_counts, status, &
                                                   relative_angular_deviations, threshold, is_outlier, gene_status, &
                                                   tmp_angular_deviations, tmp_relative_angular_deviations_perm, &
                                                   quantile_level, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), intent(in) :: signed_angles(n_genes)
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! DM_MIN(above(-PI))
            !! DM_MAX(PI)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! M_GENE_TO_FAM_DOC(signed_angles)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), intent(out) :: family_mean_angles(n_families)
            !! Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            !! `CM_FAMILY_MEAN_ANGLE_SENTINEL` where the family has fewer than three members or no
            !! stable direction (see `status`)
        real(real64), intent(out) :: angular_dispersions(n_families)
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `CM_ANGULAR_DISPERSION_SENTINEL` where `status` reports why there is none
        integer(int32), intent(out) :: member_counts(n_families)
            !! Number of genes of each family -- the members the statistics are taken over
        integer(int32), intent(out) :: status(n_families)
            !! CM_FAMILY_STATUS_DOC(mean angle, fewer than three genes belong to the family)
        real(real64), intent(out) :: relative_angular_deviations(n_genes)
            !! Absolute angular distance between each gene's angle and its family's mean angle,
            !! divided by the family's angular dispersion; `CM_RELATIVE_ANGULAR_DEVIATIONS_SENTINEL`
            !! where `gene_status` is not zero
        real(real64), intent(out) :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), intent(out) :: is_outlier(n_genes)
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(int32), intent(out) :: gene_status(n_genes)
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family, and
            !! otherwise its family's `status`
        real(real64), intent(out) :: tmp_angular_deviations(n_genes)
            !! Work array: the absolute angular distance between each gene's angle and its family's
            !! mean angle
        integer(int32), intent(out) :: tmp_relative_angular_deviations_perm(n_genes)
            !! Work array: the indices of the relative angular deviations that exist, sorted by value
        real(real64), intent(in), optional :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! DM_DEFAULT(CM_OUTLIER_QUANTILE_LEVEL_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(in), optional :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! DM_DEFAULT(CM_MIN_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! DM_DEFAULT(CM_MAX_ANGULAR_DISPERSION_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(CM_MAX_ANGULAR_DISPERSION_LIMIT)
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.

        integer(int32) :: i_gene

        call compute_family_direction_rap_impl(n_genes, n_families, signed_angles, gene_to_fam, &
                                               family_mean_angles, angular_dispersions, member_counts, status, &
                                               min_angular_dispersion, max_angular_dispersion)
        call compute_angular_deviations_rap_impl(n_genes, n_families, signed_angles, family_mean_angles, &
                                                 gene_to_fam, tmp_angular_deviations)
        call compute_relative_angular_deviations_impl(n_genes, n_families, tmp_angular_deviations, angular_dispersions, &
                                                      gene_to_fam, relative_angular_deviations)
        call compute_angle_outlier_threshold_impl(n_genes, relative_angular_deviations, &
                                                  tmp_relative_angular_deviations_perm, threshold, quantile_level)
        call flag_angle_outliers_impl(n_genes, relative_angular_deviations, threshold, is_outlier)

        ! every signed angle is a direction, so no gene lacks one here
        do concurrent (i_gene = 1:n_genes) shared(gene_status, gene_to_fam, status)
            gene_status(i_gene) = gene_status_of(gene_to_fam(i_gene), status, .true._c_bool)
        end do
    end subroutine detect_angle_outliers_rap_impl

    ! ------------------------------------------------------------------------------------------
    ! Private helpers
    ! ------------------------------------------------------------------------------------------

    !> The largest angular dispersion rounding alone can produce in a family of identical or
    !| proportional members: \((m + 1.5\,n + 16)\,\varepsilon\) for \(m\) members in \(n\)
    !| axes, with \(\varepsilon\) the machine epsilon; signed angles pass \(n = 0\).
    !|
    !| Derivation, spherical. A member's unit vector is its components times a power of two,
    !| divided by the rounded length: the sum of \(n\) squares is off by at most \(n\varepsilon\)
    !| relative, its square root by \(n\varepsilon/2 + \varepsilon\), the division adds
    !| \(2\varepsilon\) per component, so \(|u - w| \le (n/2 + 5)\varepsilon\) for the true
    !| unit vector \(w\). The family's sum of \(m\) such vectors is off by at most
    !| \((m - 1)\varepsilon\) relative per component, in any summation order, so the sum divided
    !| by \(m\) lies within \((m + n/2 + 4)\varepsilon\) of \(w\); normalising it moves it no
    !| further than that, plus the same \((n/2 + 4)\varepsilon\) of its own rounding:
    !| \(|v - w| \le (m + n + 8)\varepsilon\). The deviation \(2\,\mathrm{atan2}(|u - v|, |u + v|)\)
    !| is then at most \(|u - w| + |v - w| \le (m + 1.5\,n + 13)\varepsilon\), its own evaluation
    !| adding only relative error, and the dispersion of such deviations is their root mean square.
    !|
    !| RAP. The sums of sines and cosines are off by at most \((m - 1)\varepsilon\) relative, which
    !| turns the atan2 of them by at most as much; the library sine and cosine add \(2\varepsilon\),
    !| the atan2 an ulp of \(\pi\) (\(2\varepsilon\)), the difference to the mean and its wrapping
    !| around the circle up to \(6\varepsilon\): \((m + 9)\varepsilon\).
    !|
    !| The constant 16 covers both remainders. Measured, identical and proportional families stay
    !| below about a tenth of this in \(m\) and a hundredth in \(n\), under gfortran and under
    !| ifx's fast floating-point model alike.
    pure real(real64) function dispersion_noise_floor(n_members, n_axes) result(noise_floor)
        integer(int32), intent(in) :: n_members
            !! Number of members the dispersion is taken over
        integer(int32), intent(in) :: n_axes
            !! Number of axes of the vectors, 0 for signed angles

        noise_floor = (real(n_members, real64) + 1.5_real64*real(n_axes, real64) + 16.0_real64)*epsilon(1.0_real64)
    end function dispersion_noise_floor

    !> Angular dispersion and status of every family whose `status` is still zero, from the sum of
    !| its members' \(1 - \cos\delta\); a family already reported gets the dispersion sentinel.
    !|
    !| \(1 - R\) is that sum's mean, so \(\sigma = \sqrt{-2 \ln(1 - \text{mean})}\). A mean of 1
    !| or more is a resultant rounded to nothing, an unbounded dispersion; a dispersion within the
    !| noise floor is rounding, not variation.
    pure subroutine classify_dispersions(n_families, member_counts, n_axes, min_dispersion, max_dispersion, &
                                         angular_dispersions, status)
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        integer(int32), intent(in) :: member_counts(n_families)
            !! Number of members summed per family, at least one where `status` is zero
        integer(int32), intent(in) :: n_axes
            !! Number of axes of the vectors, 0 for signed angles
        real(real64), intent(in) :: min_dispersion
            !! Smallest accepted dispersion
        real(real64), intent(in) :: max_dispersion
            !! Largest accepted dispersion
        real(real64), intent(inout) :: angular_dispersions(n_families)
            !! On entry, per family the sum over its members of \(1 - \cos\delta\), each deviation
            !! from the family's direction; on exit the dispersion where it is accepted, otherwise
            !! `CM_ANGULAR_DISPERSION_SENTINEL`
        integer(int32), intent(inout) :: status(n_families)
            !! Zero for a family with statistics; on exit, why a family has no accepted dispersion

        real(real64) :: one_minus_resultant, sigma
        integer(int32) :: i_family

        do concurrent (i_family = 1:n_families) local(one_minus_resultant, sigma) &
            shared(member_counts, n_axes, min_dispersion, max_dispersion, angular_dispersions, status)

            if (is_err(status(i_family))) then
                angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                cycle
            end if
            one_minus_resultant = angular_dispersions(i_family)/real(member_counts(i_family), real64)
            angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
            if (one_minus_resultant >= 1.0_real64) then
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if

            sigma = sqrt(-2.0_real64*log_one_minus(one_minus_resultant))
            if (sigma <= dispersion_noise_floor(member_counts(i_family), n_axes) .or. sigma < min_dispersion) then
                call set_err(status(i_family), STAT_NO_ANGULAR_VARIATION)
            else if (sigma > max_dispersion) then
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
            else
                angular_dispersions(i_family) = sigma
            end if
        end do
    end subroutine classify_dispersions

    !> The status of one gene: why it has no relative angular deviation, or zero if it has one.
    pure integer(int32) function gene_status_of(family, family_status, gene_has_direction) result(gene_status)
        integer(int32), intent(in) :: family
            !! The gene's family, or `M_GENE_TO_FAM_SENTINEL`
        integer(int32), intent(in) :: family_status(:)
            !! The status of every family
        logical(c_bool), intent(in) :: gene_has_direction
            !! Whether the gene has a direction

        if (family == M_GENE_TO_FAM_SENTINEL) then
            gene_status = STAT_NO_FAMILY
        else if (.not. gene_has_direction) then
            gene_status = STAT_ZERO_VECTOR
        else
            gene_status = family_status(family)
        end if
    end function gene_status_of
end module tox_get_outliers_by_angle_impl
