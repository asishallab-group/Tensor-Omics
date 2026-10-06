#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_get_outliers_by_angle(module)]]
!| Gene outliers, from how far each gene's direction deviates from its family's.
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
module tox_get_outliers_by_angle_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_family_direction_c
    public :: compute_angular_deviations_c
    public :: compute_family_direction_rap_c
    public :: compute_angular_deviations_rap_c
    public :: compute_relative_angular_deviations_c
    public :: compute_angle_outlier_threshold_c
    public :: compute_angle_outlier_threshold_expert_c
    public :: flag_angle_outliers_c
    public :: detect_angle_outliers_c
    public :: detect_angle_outliers_expert_c
    public :: detect_angle_outliers_rap_c
    public :: detect_angle_outliers_rap_expert_c

contains

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_family_direction(subroutine)]]
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
    subroutine compute_family_direction_c(&
            n_axes,&
            n_genes,&
            n_families,&
            expression_vectors,&
            gene_to_fam,&
            family_directions,&
            angular_dispersions,&
            member_counts,&
            status,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="compute_family_direction_c")
        use tox_get_outliers_by_angle, only: compute_family_direction

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! Expression vector of each gene, one per column
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: family_directions
            !! Mean direction of each family as a unit vector, one per column; the zero vector
            !! where the family has fewer than three members with a direction or no stable
            !! direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family that have a direction, i.e. whose expression vector
            !! is not zero -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes of the family have a direction: no direction, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_directions, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)

        call compute_family_direction(&
            n_axes = n_axes,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            gene_to_fam = gene_to_fam,&
            family_directions = family_directions,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine compute_family_direction_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_angular_deviations(subroutine)]]
    !| The angle lies in \([0, \pi]\) and does not depend on the length of either vector.
    subroutine compute_angular_deviations_c(&
            n_axes,&
            n_genes,&
            n_families,&
            expression_vectors,&
            family_directions,&
            gene_to_fam,&
            angular_deviations,&
            ierr&
        ) bind(C, name="compute_angular_deviations_c")
        use tox_get_outliers_by_angle, only: compute_angular_deviations

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! Expression vector of each gene, one per column
        real(c_double), dimension(n_axes, n_families), intent(in), target :: family_directions
            !! Mean direction of each family, one per column, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_family_direction_impl(subroutine)]]
            !! gives it; the zero vector marks a family without one
            !! The minimum valid value is `-1.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_genes), intent(out), target :: angular_deviations
            !! Angle in radians between each gene and its family's direction, in \([0, \pi]\), or
            !! `-1.0_real64` where the gene has no family, its expression vector
            !! is zero or all subnormal, or its family has no direction
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(family_directions, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(angular_deviations, n_genes)

        call compute_angular_deviations(&
            n_axes = n_axes,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            family_directions = family_directions,&
            gene_to_fam = gene_to_fam,&
            angular_deviations = angular_deviations,&
            ierr = ierr&
        )
    end subroutine compute_angular_deviations_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_family_direction_rap(subroutine)]]
    !| Every angle is taken as a unit vector \((\cos\theta, \sin\theta)\). A family's mean angle is
    !| the angle of their sum, and with \(R\) the length of that sum divided by the number of
    !| members, its angular dispersion is \(\sigma = \sqrt{-2 \ln R}\): zero when all members have
    !| the same angle, growing as they spread. It is computed from the members' own angular
    !| deviations -- the distances `compute_angular_deviations_rap` reports -- as
    !| \(1 - R = \frac{1}{m}\sum_j 2\sin^2(\delta_j/2)\), which loses nothing to cancellation when
    !| the members are close.
    subroutine compute_family_direction_rap_c(&
            n_genes,&
            n_families,&
            signed_angles,&
            gene_to_fam,&
            family_mean_angles,&
            angular_dispersions,&
            member_counts,&
            status,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="compute_family_direction_rap_c")
        use tox_get_outliers_by_angle, only: compute_family_direction_rap

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: signed_angles
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! The minimum valid value is `above(-PI)`.
            !! The maximum valid value is `PI`.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_families), intent(out), target :: family_mean_angles
            !! Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            !! `-4.0_real64` where the family has fewer than three members or no
            !! stable direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes belong to the family: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_mean_angles, n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)

        call compute_family_direction_rap(&
            n_genes = n_genes,&
            n_families = n_families,&
            signed_angles = signed_angles,&
            gene_to_fam = gene_to_fam,&
            family_mean_angles = family_mean_angles,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine compute_family_direction_rap_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_angular_deviations_rap(subroutine)]]
    !| The distance is taken around the circle, so it lies in \([0, \pi]\): an angle just below
    !| \(\pi\) and one just above \(-\pi\) are close. See [[f42_math_impl(module):wrap_angle(function)]].
    subroutine compute_angular_deviations_rap_c(&
            n_genes,&
            n_families,&
            signed_angles,&
            family_mean_angles,&
            gene_to_fam,&
            angular_deviations,&
            ierr&
        ) bind(C, name="compute_angular_deviations_rap_c")
        use tox_get_outliers_by_angle, only: compute_angular_deviations_rap

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: signed_angles
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! The minimum valid value is `above(-PI)`.
            !! The maximum valid value is `PI`.
        real(c_double), dimension(n_families), intent(in), target :: family_mean_angles
            !! Circular mean angle of each family in radians, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_family_direction_rap_impl(subroutine)]]
            !! gives it, or `-4.0_real64` for a family without one
            !! The minimum valid value is `above(-PI)`.
            !! The maximum valid value is `PI`.
            !! The value `-4.0_real64` is additionally accepted.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_genes), intent(out), target :: angular_deviations
            !! Absolute angular distance in radians between each gene's angle and its family's mean
            !! angle, in \([0, \pi]\), or `-1.0_real64` where the gene has no
            !! family or its family has no mean angle
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_mean_angles, n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(angular_deviations, n_genes)

        call compute_angular_deviations_rap(&
            n_genes = n_genes,&
            n_families = n_families,&
            signed_angles = signed_angles,&
            family_mean_angles = family_mean_angles,&
            gene_to_fam = gene_to_fam,&
            angular_deviations = angular_deviations,&
            ierr = ierr&
        )
    end subroutine compute_angular_deviations_rap_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_relative_angular_deviations(subroutine)]]
    !| The relative angular deviation says how far off a gene is in units of its own family's
    !| spread, so genes of tight and loose families can be ranked together. It is a plain quotient,
    !| with nothing subtracted. It serves both variants: the angular deviations may come from
    !| [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_impl(subroutine)]] or
    !| from [[tox_get_outliers_by_angle_impl(module):compute_angular_deviations_rap_impl(subroutine)]].
    subroutine compute_relative_angular_deviations_c(&
            n_genes,&
            n_families,&
            angular_deviations,&
            angular_dispersions,&
            gene_to_fam,&
            relative_angular_deviations,&
            ierr&
        ) bind(C, name="compute_relative_angular_deviations_c")
        use tox_get_outliers_by_angle, only: compute_relative_angular_deviations

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: angular_deviations
            !! Angular deviation of each gene from its family in radians, or
            !! `-1.0_real64` where there is none
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `PI`.
            !! The value `-1.0_real64` is additionally accepted.
        real(c_double), dimension(n_families), intent(in), target :: angular_dispersions
            !! Angular dispersion of each family, or `-1.0_real64` where there
            !! is none
            !! The minimum valid value is `0.0_real64`.
            !! The value `-1.0_real64` is additionally accepted.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `angular_deviations`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_genes), intent(out), target :: relative_angular_deviations
            !! Angular deviation divided by the family's angular dispersion, at least zero, or
            !! `-1.0_real64` where the gene has no family or no
            !! angular deviation, or its family's dispersion is the sentinel, zero, or below the
            !! smallest normal number `tiny(1.0)` (about `2.2e-308`) -- a dispersion so small that
            !! the quotient could exceed the largest representable number
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_ARRAY_NON_NULL(angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)

        call compute_relative_angular_deviations(&
            n_genes = n_genes,&
            n_families = n_families,&
            angular_deviations = angular_deviations,&
            angular_dispersions = angular_dispersions,&
            gene_to_fam = gene_to_fam,&
            relative_angular_deviations = relative_angular_deviations,&
            ierr = ierr&
        )
    end subroutine compute_relative_angular_deviations_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_angle_outlier_threshold(subroutine)]]
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
    subroutine compute_angle_outlier_threshold_c(&
            n_genes,&
            relative_angular_deviations,&
            threshold,&
            quantile_level,&
            ierr&
        ) bind(C, name="compute_angle_outlier_threshold_c")
        use tox_get_outliers_by_angle, only: compute_angle_outlier_threshold

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        real(c_double), dimension(n_genes), intent(in), target :: relative_angular_deviations
            !! Relative angular deviation of each gene, or `-1.0_real64`
            !! where there is none
            !! The minimum valid value is `0.0_real64`.
            !! The value `-1.0_real64` is additionally accepted.
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)

        call compute_angle_outlier_threshold(&
            n_genes = n_genes,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            quantile_level = quantile_level,&
            ierr = ierr&
        )
    end subroutine compute_angle_outlier_threshold_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):compute_angle_outlier_threshold_expert(subroutine)]]
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
    subroutine compute_angle_outlier_threshold_expert_c(&
            n_genes,&
            relative_angular_deviations,&
            tmp_relative_angular_deviations_perm,&
            threshold,&
            quantile_level,&
            ierr&
        ) bind(C, name="compute_angle_outlier_threshold_expert_c")
        use tox_get_outliers_by_angle, only: compute_angle_outlier_threshold_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        real(c_double), dimension(n_genes), intent(in), target :: relative_angular_deviations
            !! Relative angular deviation of each gene, or `-1.0_real64`
            !! where there is none
            !! The minimum valid value is `0.0_real64`.
            !! The value `-1.0_real64` is additionally accepted.
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_relative_angular_deviations_perm
            !! Work array: the indices of the values that exist, sorted by value
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_relative_angular_deviations_perm, n_genes)

        call compute_angle_outlier_threshold_expert(&
            n_genes = n_genes,&
            relative_angular_deviations = relative_angular_deviations,&
            tmp_relative_angular_deviations_perm = tmp_relative_angular_deviations_perm,&
            threshold = threshold,&
            quantile_level = quantile_level,&
            ierr = ierr&
        )
    end subroutine compute_angle_outlier_threshold_expert_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):flag_angle_outliers(subroutine)]]
    !| A gene is an outlier when its relative angular deviation is positive and strictly greater
    !| than `threshold`; a value equal to the threshold is not flagged. So a gene without a value,
    !| or with a deviation of zero, is never flagged, whatever the threshold.
    subroutine flag_angle_outliers_c(&
            n_genes,&
            relative_angular_deviations,&
            threshold,&
            is_outlier,&
            ierr&
        ) bind(C, name="flag_angle_outliers_c")
        use tox_get_outliers_by_angle, only: flag_angle_outliers

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        real(c_double), dimension(n_genes), intent(in), target :: relative_angular_deviations
            !! Relative angular deviation of each gene, or `-1.0_real64`
            !! where there is none
            !! The minimum valid value is `0.0_real64`.
            !! The value `-1.0_real64` is additionally accepted.
        real(c_double), intent(in), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier, as
            !! [[tox_get_outliers_by_angle_impl(module):compute_angle_outlier_threshold_impl(subroutine)]]
            !! computes it
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! `.true.` for each gene that is an outlier
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)

        call flag_angle_outliers(&
            n_genes = n_genes,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            is_outlier = is_outlier,&
            ierr = ierr&
        )
    end subroutine flag_angle_outliers_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):detect_angle_outliers(subroutine)]]
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
    subroutine detect_angle_outliers_c(&
            n_axes,&
            n_genes,&
            n_families,&
            expression_vectors,&
            gene_to_fam,&
            family_directions,&
            angular_dispersions,&
            member_counts,&
            status,&
            relative_angular_deviations,&
            threshold,&
            is_outlier,&
            gene_status,&
            quantile_level,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="detect_angle_outliers_c")
        use tox_get_outliers_by_angle, only: detect_angle_outliers

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! Expression vector of each gene, one per column
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: family_directions
            !! Mean direction of each family as a unit vector, one per column; the zero vector
            !! where the family has fewer than three members with a direction or no stable
            !! direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family that have a direction, i.e. whose expression vector
            !! is not zero -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes of the family have a direction: no direction, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.
        real(c_double), dimension(n_genes), intent(out), target :: relative_angular_deviations
            !! Angle between each gene and its family's direction, divided by the family's angular
            !! dispersion; `-1.0_real64` where `gene_status` is not zero
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(c_int), dimension(n_genes), intent(out), target :: gene_status
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family,
            !! [[tox_errors(module):STAT_ZERO_VECTOR(variable)]] for one whose expression vector is
            !! zero or all subnormal, and otherwise its family's `status`
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_directions, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_status, n_genes)

        call detect_angle_outliers(&
            n_axes = n_axes,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            gene_to_fam = gene_to_fam,&
            family_directions = family_directions,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            is_outlier = is_outlier,&
            gene_status = gene_status,&
            quantile_level = quantile_level,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine detect_angle_outliers_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):detect_angle_outliers_expert(subroutine)]]
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
    subroutine detect_angle_outliers_expert_c(&
            n_axes,&
            n_genes,&
            n_families,&
            expression_vectors,&
            gene_to_fam,&
            family_directions,&
            angular_dispersions,&
            member_counts,&
            status,&
            relative_angular_deviations,&
            threshold,&
            is_outlier,&
            gene_status,&
            tmp_angular_deviations,&
            tmp_relative_angular_deviations_perm,&
            quantile_level,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="detect_angle_outliers_expert_c")
        use tox_get_outliers_by_angle, only: detect_angle_outliers_expert

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes, i.e. the length of each expression vector
        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_axes, n_genes), intent(in), target :: expression_vectors
            !! Expression vector of each gene, one per column
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_axes, n_families), intent(out), target :: family_directions
            !! Mean direction of each family as a unit vector, one per column; the zero vector
            !! where the family has fewer than three members with a direction or no stable
            !! direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family that have a direction, i.e. whose expression vector
            !! is not zero -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a direction and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes of the family have a direction: no direction, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no direction, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the direction kept.
        real(c_double), dimension(n_genes), intent(out), target :: relative_angular_deviations
            !! Angle between each gene and its family's direction, divided by the family's angular
            !! dispersion; `-1.0_real64` where `gene_status` is not zero
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(c_int), dimension(n_genes), intent(out), target :: gene_status
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family,
            !! [[tox_errors(module):STAT_ZERO_VECTOR(variable)]] for one whose expression vector is
            !! zero or all subnormal, and otherwise its family's `status`
        real(c_double), dimension(n_genes), intent(out), target :: tmp_angular_deviations
            !! Work array: the angle between each gene and its family's direction
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_relative_angular_deviations_perm
            !! Work array: the indices of the relative angular deviations that exist, sorted by value
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_directions, n_axes * n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_status, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_relative_angular_deviations_perm, n_genes)

        call detect_angle_outliers_expert(&
            n_axes = n_axes,&
            n_genes = n_genes,&
            n_families = n_families,&
            expression_vectors = expression_vectors,&
            gene_to_fam = gene_to_fam,&
            family_directions = family_directions,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            is_outlier = is_outlier,&
            gene_status = gene_status,&
            tmp_angular_deviations = tmp_angular_deviations,&
            tmp_relative_angular_deviations_perm = tmp_relative_angular_deviations_perm,&
            quantile_level = quantile_level,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine detect_angle_outliers_expert_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):detect_angle_outliers_rap(subroutine)]]
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
    subroutine detect_angle_outliers_rap_c(&
            n_genes,&
            n_families,&
            signed_angles,&
            gene_to_fam,&
            family_mean_angles,&
            angular_dispersions,&
            member_counts,&
            status,&
            relative_angular_deviations,&
            threshold,&
            is_outlier,&
            gene_status,&
            quantile_level,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="detect_angle_outliers_rap_c")
        use tox_get_outliers_by_angle, only: detect_angle_outliers_rap

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: signed_angles
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! The minimum valid value is `above(-PI)`.
            !! The maximum valid value is `PI`.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_families), intent(out), target :: family_mean_angles
            !! Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            !! `-4.0_real64` where the family has fewer than three members or no
            !! stable direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes belong to the family: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.
        real(c_double), dimension(n_genes), intent(out), target :: relative_angular_deviations
            !! Absolute angular distance between each gene's angle and its family's mean angle,
            !! divided by the family's angular dispersion; `-1.0_real64`
            !! where `gene_status` is not zero
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(c_int), dimension(n_genes), intent(out), target :: gene_status
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family, and
            !! otherwise its family's `status`
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_mean_angles, n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_status, n_genes)

        call detect_angle_outliers_rap(&
            n_genes = n_genes,&
            n_families = n_families,&
            signed_angles = signed_angles,&
            gene_to_fam = gene_to_fam,&
            family_mean_angles = family_mean_angles,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            is_outlier = is_outlier,&
            gene_status = gene_status,&
            quantile_level = quantile_level,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine detect_angle_outliers_rap_c

    !> summary: C-wrapper for [[tox_get_outliers_by_angle(module):detect_angle_outliers_rap_expert(subroutine)]]
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
    subroutine detect_angle_outliers_rap_expert_c(&
            n_genes,&
            n_families,&
            signed_angles,&
            gene_to_fam,&
            family_mean_angles,&
            angular_dispersions,&
            member_counts,&
            status,&
            relative_angular_deviations,&
            threshold,&
            is_outlier,&
            gene_status,&
            tmp_angular_deviations,&
            tmp_relative_angular_deviations_perm,&
            quantile_level,&
            min_angular_dispersion,&
            max_angular_dispersion,&
            ierr&
        ) bind(C, name="detect_angle_outliers_rap_expert_c")
        use tox_get_outliers_by_angle, only: detect_angle_outliers_rap_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes
        integer(c_int), intent(in), target :: n_families
            !! Number of gene families
        real(c_double), dimension(n_genes), intent(in), target :: signed_angles
            !! Signed angle of each gene in radians, in \((-\pi, \pi]\), as
            !! [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]]
            !! measures it. A gene without an angle is one with no family in `gene_to_fam`; its
            !! entry here is ignored, but must still lie in the range.
            !! The minimum valid value is `above(-PI)`.
            !! The maximum valid value is `PI`.
        integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `signed_angles`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(c_double), dimension(n_families), intent(out), target :: family_mean_angles
            !! Circular mean angle of each family in radians, in \((-\pi, \pi]\), or
            !! `-4.0_real64` where the family has fewer than three members or no
            !! stable direction (see `status`)
        real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
            !! Angular dispersion \(\sigma = \sqrt{-2 \ln R}\) of each family, or
            !! `-1.0_real64` where `status` reports why there is none
        integer(c_int), dimension(n_families), intent(out), target :: member_counts
            !! Number of genes of each family -- the members the statistics are taken over
        integer(c_int), dimension(n_families), intent(out), target :: status
            !! Zero where the family has a mean angle and an accepted angular dispersion, otherwise the reason it has not. [[tox_errors(module):STAT_TOO_FEW_MEMBERS(variable)]] where fewer than three genes belong to the family: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_STABLE_DIRECTION(variable)]] where the unit vectors cancel out, exactly or up to rounding, or the dispersion exceeds `max_angular_dispersion`: no mean angle, no dispersion. [[tox_errors(module):STAT_NO_ANGULAR_VARIATION(variable)]] where the dispersion is too small to tell from rounding, whatever `min_angular_dispersion` says, or below `min_angular_dispersion`, with the mean angle kept.
        real(c_double), dimension(n_genes), intent(out), target :: relative_angular_deviations
            !! Absolute angular distance between each gene's angle and its family's mean angle,
            !! divided by the family's angular dispersion; `-1.0_real64`
            !! where `gene_status` is not zero
        real(c_double), intent(out), target :: threshold
            !! The relative angular deviation a gene must exceed to be an outlier; the largest
            !! representable number when fewer than two genes have a relative angular deviation
        logical(c_bool), dimension(n_genes), intent(out), target :: is_outlier
            !! `.true.` for each gene whose relative angular deviation is positive and strictly above
            !! `threshold`
        integer(c_int), dimension(n_genes), intent(out), target :: gene_status
            !! Zero where the gene has a relative angular deviation, otherwise the reason it has
            !! not: [[tox_errors(module):STAT_NO_FAMILY(variable)]] for a gene without a family, and
            !! otherwise its family's `status`
        real(c_double), dimension(n_genes), intent(out), target :: tmp_angular_deviations
            !! Work array: the absolute angular distance between each gene's angle and its family's
            !! mean angle
        integer(c_int), dimension(n_genes), intent(out), target :: tmp_relative_angular_deviations_perm
            !! Work array: the indices of the relative angular deviations that exist, sorted by value
        real(c_double), intent(in), target :: quantile_level
            !! Quantile level in \([0, 1]\) the threshold is taken at
            !! The default value is `0.95_real64`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: min_angular_dispersion
            !! Smallest accepted angular dispersion. Where it exceeds `max_angular_dispersion`,
            !! no dispersion is accepted and every family with statistics is reported.
            !! The default value is `0.0_real64`.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: max_angular_dispersion
            !! Largest accepted angular dispersion; the default is the dispersion of \(R = 1/2\)
            !! The default value is `sqrt(-2.0_real64*log(0.5_real64))`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `5.0_real64`.
            !! At most 5 (\(R \ge e^{-12.5} \approx 4 \times 10^{-6}\)): above that, a resultant that cancels only up to rounding would pass as a stable direction.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_families)
        M_CHECK_NON_NULL(threshold)
        M_CHECK_NON_NULL(quantile_level)
        M_CHECK_NON_NULL(min_angular_dispersion)
        M_CHECK_NON_NULL(max_angular_dispersion)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_to_fam, n_genes)
        M_CHECK_ARRAY_NON_NULL(family_mean_angles, n_families)
        M_CHECK_ARRAY_NON_NULL(angular_dispersions, n_families)
        M_CHECK_ARRAY_NON_NULL(member_counts, n_families)
        M_CHECK_ARRAY_NON_NULL(status, n_families)
        M_CHECK_ARRAY_NON_NULL(relative_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(is_outlier, n_genes)
        M_CHECK_ARRAY_NON_NULL(gene_status, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_angular_deviations, n_genes)
        M_CHECK_ARRAY_NON_NULL(tmp_relative_angular_deviations_perm, n_genes)

        call detect_angle_outliers_rap_expert(&
            n_genes = n_genes,&
            n_families = n_families,&
            signed_angles = signed_angles,&
            gene_to_fam = gene_to_fam,&
            family_mean_angles = family_mean_angles,&
            angular_dispersions = angular_dispersions,&
            member_counts = member_counts,&
            status = status,&
            relative_angular_deviations = relative_angular_deviations,&
            threshold = threshold,&
            is_outlier = is_outlier,&
            gene_status = gene_status,&
            tmp_angular_deviations = tmp_angular_deviations,&
            tmp_relative_angular_deviations_perm = tmp_relative_angular_deviations_perm,&
            quantile_level = quantile_level,&
            min_angular_dispersion = min_angular_dispersion,&
            max_angular_dispersion = max_angular_dispersion,&
            ierr = ierr&
        )
    end subroutine detect_angle_outliers_rap_expert_c

end module tox_get_outliers_by_angle_c
#endif
