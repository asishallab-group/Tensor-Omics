#include <src/macros.h>

!> Normalized tissue (axis) versatility: how uniformly a gene is expressed across tissues.
!|
!| An angle-based metric. A gene expressed equally across every selected axis points along the
!| diagonal of that subspace and scores maximally versatile; one confined to a single tissue
!| points along that axis and scores minimally. Normalized, so scores over different numbers of
!| selected axes are comparable.
module tox_tissue_versatility_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_math_impl, only: degrees
    M_IMPLICIT_NONE
contains

    !> summary: Computes normalized tissue versatility for selected expression vectors.
    !| AUTHOR_VIVIAN_BASS
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    !|
    !| The masks follow the `n_selected_` convention, so the generated wrapper validates that each
    !| selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.
    pure subroutine compute_tissue_versatility_impl(n_axes, n_vectors, expression_vectors, vectors_selection_mask, &
                                               n_selected_vectors, axes_selection_mask, n_selected_axes, &
                                               tissue_versatilities, tissue_angles_deg)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(int32), intent(in) :: n_vectors
            !! Number of expression vectors (genes)
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection_mask)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in vectors_selection_mask)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical(c_bool), dimension(n_vectors), intent(in) :: vectors_selection_mask
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical(c_bool), dimension(n_axes), intent(in) :: axes_selection_mask
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_versatilities
            !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_angles_deg
            !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees

        ! Local variables
        integer(int32) :: i_vec, i_axis, out_idx, first_axis, scale_exponent
        real(real64) :: norm_diag, norm_factor, max_abs, scale_factor, pivot, scaled_value, mean_offset_from_pivot
        real(real64) :: scaled_sum, scaled_norm_sq, rejection_sq, projection, scaled_norm, cos_phi, one_minus_cos

        ! Handle edge case: when only one axis is selected, tissue versatility is always 0
        if (n_selected_axes == 1) then
            tissue_versatilities = 0.0_real64
            tissue_angles_deg = 0.0_real64
            return
        end if

        norm_diag = sqrt(real(n_selected_axes, real64))
        ! cos_phi against the space diagonal ranges from 1 (uniform expression across all selected
        ! axes) down to 1/norm_diag (all expression concentrated in a single axis), so (1 - cos_phi)
        ! ranges over [0, 1 - 1/norm_diag]. Dividing by that upper bound rescales versatility to [0, 1].
        norm_factor = 1.0_real64 - 1.0_real64/norm_diag
        ! Not findloc(axes_selection_mask, .true._c_bool): ifx stores .true. as all bits set, while
        ! C, numpy and R pass the byte 1 in a logical(c_bool). findloc, like .eqv., compares the
        ! whole byte, so it missed a true from C, returned 0, and the read out of bounds segfaulted
        ! under ifx at every optimization level. Testing the mask itself works on both compilers:
        ! ifx looks at the low bit, gfortran at any nonzero byte.
        first_axis = 1
        do while (.not. axes_selection_mask(first_axis))
            first_axis = first_axis + 1
        end do

        ! Each gene v is split into its projection on the unit space diagonal, p = sum(v)/norm_diag,
        ! and its rejection from it, whose length r is that of v minus its mean. The angle is
        ! atan2(r, p), which unlike acos(cos_phi) keeps its precision near 0.
        ! 1 - cos_phi comes from one of two formulas, split at cos_phi = 3/4:
        ! - above 3/4, r**2/(|v|*(|v| + p)), which avoids the cancellation of 1 - cos_phi near 1;
        ! - at or below 3/4, 1 - cos_phi, which is then at least 1/4, so an error in cos_phi grows
        !   by at most cos_phi/(1 - cos_phi) = 3.
        ! The split is there for the single-axis gene, not for accuracy: both formulas are about
        ! equally accurate between 1/2 and 3/4 (around 10 ulps of TV against real128), but only
        ! the second makes a single-axis gene exact. Any split above 1/sqrt(2) would do.
        ! cos_phi is computed as (sum(v)/|v|)/norm_diag. For a single-axis gene sum(v)/|v| is
        ! exactly 1, since |v| = sqrt(v_k**2) is |v_k| exactly, so cos_phi is 1/norm_diag computed
        ! as in norm_factor, and TV is exactly 1 over any number of axes: its cos_phi, 1/sqrt(n),
        ! is at most 1/sqrt(2), below the split. r**2 is a rounded sum, and 1 - p/|v| rounds p
        ! first; either left many single-axis values an ulp or two off 1, some above it.
        ! The rejection is taken of v minus its first selected value (a multiple of the diagonal,
        ! so the rejection is the same), which is exactly zero for a uniform gene: angle and TV
        ! are then exactly 0.
        ! Every value is multiplied by scale_factor = 2**-scale_exponent, which is exact unless the
        ! product underflows. 2**scale_exponent is the power of two just above the gene's largest
        ! magnitude, so the largest scaled value lies in [1/2, 1): the sums of squares cannot
        ! overflow, and the largest square cannot underflow. A component far below the largest can
        ! still underflow when scaled or squared; it then adds less to its sum than that sum's
        ! rounding. scale_exponent is clamped at both ends, which keeps scale_factor a normal number:
        ! - at 1 - maxexponent = -1023 from below: for a largest magnitude below 2**-1024 (deep in
        !   the subnormal range, which starts below 2**-1022), 2**-scale_exponent would overflow.
        !   Such a gene is scaled up by 2**1023, exactly, and its largest value lands in
        !   [2**-51, 1/2), still far from both limits.
        ! - at maxexponent - 2 = 1022 from above: for a largest magnitude of 2**1022 or more,
        !   2**-scale_exponent would be subnormal, which flush-to-zero (a caller built with fast
        !   math may leave it on) turns into 0, and TV into NaN. Such a gene is scaled by 2**-1022, and
        !   its largest value lands in [1, 4): every square and every term of r**2 is below 100, so
        !   no sum overflows for any number of axes an int32 can count.
        ! Scaling by a power of two changes no rounding, except which tiny components underflow,
        ! so the gene gets the same result as if it were scaled into [1/2, 1).
        ! A gene that is zero in all selected axes has no direction; it gets TV 1 (maximum
        ! specificity) and an angle of 90 degrees.
        !TODO optimize: this per-gene loop is data-parallel (each i_vec is independent) but is forced sequential here only
        ! because out_idx is an accumulating compaction counter. Since n_selected_vectors is already known, this
        ! could be parallelized with `do concurrent` by precomputing an exclusive prefix sum of
        ! vectors_selection_mask to get each vector's output slot directly, consistent with how other modules in
        ! this codebase parallelize per-gene work.
        out_idx = 0
        do i_vec = 1, n_vectors
            if (.not. vectors_selection_mask(i_vec)) cycle  ! Skip if not selected
            out_idx = out_idx + 1

            max_abs = 0.0_real64
            do i_axis = 1, n_axes
                if (axes_selection_mask(i_axis)) max_abs = max(max_abs, abs(expression_vectors(i_axis, i_vec)))
            end do
            if (max_abs <= 0.0_real64) then
                ! exactly zero, not "small": a gene expressed at a tiny scale still has a direction
                tissue_versatilities(out_idx) = 1.0_real64
                tissue_angles_deg(out_idx) = 90.0_real64
                cycle
            end if
            scale_exponent = min(max(exponent(max_abs), 1 - maxexponent(max_abs)), maxexponent(max_abs) - 2)
            scale_factor = scale(1.0_real64, -scale_exponent)
            pivot = expression_vectors(first_axis, i_vec)*scale_factor

            scaled_sum = 0.0_real64
            scaled_norm_sq = 0.0_real64
            mean_offset_from_pivot = 0.0_real64
            do i_axis = 1, n_axes
                if (.not. axes_selection_mask(i_axis)) cycle
                scaled_value = expression_vectors(i_axis, i_vec)*scale_factor
                scaled_sum = scaled_sum + scaled_value
                scaled_norm_sq = scaled_norm_sq + scaled_value**2
                mean_offset_from_pivot = mean_offset_from_pivot + (scaled_value - pivot)
            end do
            mean_offset_from_pivot = mean_offset_from_pivot/real(n_selected_axes, real64)

            rejection_sq = 0.0_real64
            do i_axis = 1, n_axes
                if (.not. axes_selection_mask(i_axis)) cycle
                scaled_value = expression_vectors(i_axis, i_vec)*scale_factor
                rejection_sq = rejection_sq + ((scaled_value - pivot) - mean_offset_from_pivot)**2
            end do

            projection = scaled_sum/norm_diag
            scaled_norm = sqrt(scaled_norm_sq)
            cos_phi = (scaled_sum/scaled_norm)/norm_diag
            if (cos_phi > 0.75_real64) then
                one_minus_cos = rejection_sq/(scaled_norm*(scaled_norm + projection))
            else
                one_minus_cos = 1.0_real64 - cos_phi
            end if
            tissue_versatilities(out_idx) = one_minus_cos/norm_factor
            tissue_angles_deg(out_idx) = degrees(atan2(sqrt(rejection_sq), projection))
        end do

    end subroutine compute_tissue_versatility_impl

end module tox_tissue_versatility_impl
