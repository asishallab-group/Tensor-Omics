!> Unit test suite for tox_get_outliers_by_angle, the angle-based outlier detection, and for
!| the `wrap_angle` it builds on.
!|
!| Expected values are derived by hand in the comments, from inputs chosen so the angles are
!| exact: axis-aligned, orthogonal and antipodal vectors, and `cos = 1/2` for 60 degrees. A
!| tolerance appears only where a transcendental function (`acos`, `log`, `sqrt`, `atan2`, `sin`)
!| or the rounded constant `PI` stands between the input and the value, and is then a few units
!| in the last place of an O(1) result. Outputs are pre-filled with `UNWRITTEN` values so an
!| entry the procedure forgot to write cannot pass by accident.
module mod_test_tox_get_outliers_by_angle
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_support_underflow_control, &
                                             ieee_get_underflow_mode, ieee_set_underflow_mode
    use f42_math_impl, only: PI, wrap_angle
    use tox_get_outliers_by_angle
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_tox_get_outliers_by_angle

    !> A few units in the last place of an O(1) value: the rounding of one transcendental
    !| function call, or of `PI` itself
    real(real64), parameter :: TOL_ULP = 1.0e-15_real64
    !> Where several rounded operations stack up (a sum of sines, a logarithm of a rounded sum)
    real(real64), parameter :: TOL_CHAIN = 1.0e-13_real64

    !> A quiet NaN and positive infinity by their IEEE 754 bit patterns, as constants: an array
    !| constructor holding them then needs no temporary
    real(real64), parameter :: NAN_VALUE = transfer(int(z'7FF8000000000000', int64), 1.0_real64)
    real(real64), parameter :: INF_VALUE = transfer(int(z'7FF0000000000000', int64), 1.0_real64)

    real(real64), parameter :: UNWRITTEN = 777.0_real64
    integer(int32), parameter :: UNWRITTEN_INT = -777_int32

    ! The sentinels and defaults the implementation documents, repeated here on purpose: they are
    ! the published contract, and a change to one in the implementation must fail this suite
    real(real64), parameter :: ANGULAR_DEVIATIONS_SENTINEL = -1.0_real64
    real(real64), parameter :: ANGULAR_DISPERSION_SENTINEL = -1.0_real64
    real(real64), parameter :: RELATIVE_SENTINEL = -1.0_real64
    real(real64), parameter :: FAMILY_MEAN_ANGLE_SENTINEL = -4.0_real64

    real(real64), parameter :: MIN_DISPERSION_DEFAULT = 0.0_real64
    real(real64), parameter :: MAX_DISPERSION_DEFAULT = sqrt(-2.0_real64*log(0.5_real64))
    real(real64), parameter :: MAX_DISPERSION_LIMIT = 5.0_real64

contains

    function get_all_tests_tox_get_outliers_by_angle() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(41))
        all_tests(1) = test_case("test_wrap_angle", test_wrap_angle)
        all_tests(2) = test_case("test_family_direction_exact", test_family_direction_exact)
        all_tests(3) = test_case("test_family_direction_zero_gene_skipped", test_family_direction_zero_gene_skipped)
        all_tests(4) = test_case("test_family_direction_too_few_members", test_family_direction_too_few_members)
        all_tests(5) = test_case("test_family_direction_identical_vectors_never_nan", &
                                 test_family_direction_identical_vectors_never_nan)
        all_tests(6) = test_case("test_family_direction_dispersion_bounds", test_family_direction_dispersion_bounds)
        all_tests(7) = test_case("test_family_direction_extreme_magnitudes", test_family_direction_extreme_magnitudes)
        all_tests(8) = test_case("test_angular_deviations_exact", test_angular_deviations_exact)
        all_tests(9) = test_case("test_family_direction_rap_exact", test_family_direction_rap_exact)
        all_tests(10) = test_case("test_family_direction_rap_too_few_and_clamp", test_family_direction_rap_too_few_and_clamp)
        all_tests(11) = test_case("test_angular_deviations_rap_wrapping", test_angular_deviations_rap_wrapping)
        all_tests(12) = test_case("test_relative_angular_deviations_exact", test_relative_angular_deviations_exact)
        all_tests(13) = test_case("test_threshold_quantile", test_threshold_quantile)
        all_tests(14) = test_case("test_threshold_fewer_than_two_values", test_threshold_fewer_than_two_values)
        all_tests(15) = test_case("test_threshold_and_flags_with_ties_and_zeros", test_threshold_and_flags_with_ties_and_zeros)
        all_tests(16) = test_case("test_flag_rule", test_flag_rule)
        all_tests(17) = test_case("test_largest_flagged_iff_above_threshold", test_largest_flagged_iff_above_threshold)
        all_tests(18) = test_case("test_detect_angle_outliers", test_detect_angle_outliers)
        all_tests(19) = test_case("test_detect_angle_outliers_no_family_has_statistics", &
                                  test_detect_angle_outliers_no_family_has_statistics)
        all_tests(20) = test_case("test_detect_angle_outliers_rap", test_detect_angle_outliers_rap)
        all_tests(21) = test_case("test_validation_compute_family_direction", test_validation_compute_family_direction)
        all_tests(22) = test_case("test_validation_compute_angular_deviations", test_validation_compute_angular_deviations)
        all_tests(23) = test_case("test_validation_compute_family_direction_rap", test_validation_compute_family_direction_rap)
        all_tests(24) = test_case("test_validation_compute_angular_deviations_rap", &
                                  test_validation_compute_angular_deviations_rap)
        all_tests(25) = test_case("test_validation_compute_relative_angular_deviations", &
                                  test_validation_compute_relative_angular_deviations)
        all_tests(26) = test_case("test_validation_compute_angle_outlier_threshold", &
                                  test_validation_compute_angle_outlier_threshold)
        all_tests(27) = test_case("test_validation_flag_angle_outliers", test_validation_flag_angle_outliers)
        all_tests(28) = test_case("test_validation_detect_angle_outliers", test_validation_detect_angle_outliers)
        all_tests(29) = test_case("test_validation_detect_angle_outliers_rap", test_validation_detect_angle_outliers_rap)
        all_tests(30) = test_case("test_status_codes_are_distinct_statuses", test_status_codes_are_distinct_statuses)
        all_tests(31) = test_case("test_rap_exact_cancellation", test_rap_exact_cancellation)
        all_tests(32) = test_case("test_noise_floor_large_families", test_noise_floor_large_families)
        all_tests(33) = test_case("test_rap_pi_mixed_with_minus_pi", test_rap_pi_mixed_with_minus_pi)
        all_tests(34) = test_case("test_tiny_dispersion_above_and_below_floor", test_tiny_dispersion_above_and_below_floor)
        all_tests(35) = test_case("test_single_deviating_member", test_single_deviating_member)
        all_tests(36) = test_case("test_detect_identical_family_default_min", test_detect_identical_family_default_min)
        all_tests(37) = test_case("test_detect_999_plus_1_not_masked", test_detect_999_plus_1_not_masked)
        all_tests(38) = test_case("test_published_defaults", test_published_defaults)
        all_tests(39) = test_case("test_nearly_cancelling_resultant", test_nearly_cancelling_resultant)
        all_tests(40) = test_case("test_subnormal_vector_is_zero_vector", test_subnormal_vector_is_zero_vector)
        all_tests(41) = test_case("test_cancellation_never_stable_up_to_the_limit", &
                                  test_cancellation_never_stable_up_to_the_limit)
    end function get_all_tests_tox_get_outliers_by_angle

    ! ------------------------------------------------------------------------------------------
    ! helpers
    ! ------------------------------------------------------------------------------------------

    !> Assert that no entry of a real array is NaN or infinite.
    subroutine assert_all_finite(values, n, msg)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: values(n)
        character(*), intent(in) :: msg

        call assert_true(all(ieee_is_finite(values)), msg//": all finite")
    end subroutine assert_all_finite

    ! ------------------------------------------------------------------------------------------
    ! wrap_angle
    ! ------------------------------------------------------------------------------------------

    !> wrap_angle maps into (-pi, pi]: pi and -pi both give pi, and a full turn is removed.
    subroutine test_wrap_angle()
        real(real64) :: just_above_pi, wrapped

        ! PI - modulo(PI - x, 2*PI) with x = 0: modulo(PI, 2*PI) = PI, so 0 exactly
        call assert_equal_real(wrap_angle(0.0_real64), 0.0_real64, 0.0_real64, "wrap(0)")
        ! x = PI: modulo(0, 2*PI) = 0, so PI exactly
        call assert_equal_real(wrap_angle(PI), PI, 0.0_real64, "wrap(pi) stays pi")
        ! x = -PI: modulo(2*PI, 2*PI) = 0, so PI exactly -- the interval is open at -pi
        call assert_equal_real(wrap_angle(-PI), PI, 0.0_real64, "wrap(-pi) is pi")
        ! x = PI/2: modulo(PI/2, 2*PI) = PI/2, so PI/2 exactly (halving is exact)
        call assert_equal_real(wrap_angle(PI/2), PI/2, 0.0_real64, "wrap(pi/2)")
        call assert_equal_real(wrap_angle(-PI/2), -PI/2, TOL_ULP, "wrap(-pi/2)")
        ! a full turn is removed: 3*pi/2 is -pi/2, 2*pi is 0, 3*pi is pi, -3*pi/2 is pi/2
        call assert_equal_real(wrap_angle(3*PI/2), -PI/2, TOL_ULP, "wrap(3pi/2)")
        call assert_equal_real(wrap_angle(2*PI), 0.0_real64, TOL_ULP, "wrap(2pi)")
        call assert_equal_real(wrap_angle(3*PI), PI, TOL_ULP, "wrap(3pi)")
        call assert_equal_real(wrap_angle(-3*PI/2), PI/2, TOL_ULP, "wrap(-3pi/2)")

        ! an angle already in (-pi, pi] comes back unchanged, to the last bit
        call assert_equal_real(wrap_angle(-1.0e-20_real64), -1.0e-20_real64, 0.0_real64, "wrap(-1e-20) unchanged")
        call assert_equal_real(wrap_angle(nearest(-PI, 1.0_real64)), nearest(-PI, 1.0_real64), 0.0_real64, &
                               "wrap(-pi + ulp) unchanged")
        call assert_equal_real(wrap_angle(0.7_real64), 0.7_real64, 0.0_real64, "wrap(0.7) unchanged")

        ! one ulp beyond pi is the same direction as pi: rounding may land it on exactly -pi,
        ! which must come back as pi, never below or at -pi
        just_above_pi = nearest(PI, 1.0_real64)
        wrapped = wrap_angle(just_above_pi)
        call assert_true(wrapped > -PI, "wrap(pi + ulp) is inside (-pi, pi]")
        call assert_true(wrapped <= PI, "wrap(pi + ulp) is at most pi")
        call assert_equal_real(abs(wrapped), PI, TOL_ULP, "wrap(pi + ulp) is next to +-pi")
    end subroutine test_wrap_angle

    ! ------------------------------------------------------------------------------------------
    ! compute_family_direction
    ! ------------------------------------------------------------------------------------------

    !> Three families with exact statistics, and one unassigned gene.
    subroutine test_family_direction_exact()
        integer(int32), parameter :: n_axes = 2, n_genes = 12, n_families = 3
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families)
        integer(int32) :: member_counts(n_families), status(n_families), ierr

        ! family 1: (1,0), (3,0), (0,2), (0,0.5) have the unit vectors e1, e1, e2, e2. Their sum is
        !   (2,2), of length 2*sqrt(2), so R = 2*sqrt(2)/4 = sqrt(2)/2, the direction is
        !   (1,1)/sqrt(2), and sigma = sqrt(-2 ln(sqrt(2)/2)) = sqrt(ln 2) = 0.8326, accepted.
        ! family 2: (1,0), (0,1), (-1,0) sum to (0,1) exactly: R = 1/3 and
        !   sigma = sqrt(2 ln 3) = 1.4823 > 1.1774, so no stable direction: the zero vector.
        ! family 3: (1,0), (-1,0), (2,0), (-2,0) sum to exactly zero: no direction at all.
        ! gene 12 belongs to no family and must not reach any sum.
        expression_vectors = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, &
                                      0.0_real64, 2.0_real64, 0.0_real64, 0.5_real64, &
                                      1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64, &
                                      1.0_real64, 0.0_real64, -1.0_real64, 0.0_real64, &
                                      2.0_real64, 0.0_real64, -2.0_real64, 0.0_real64, &
                                      5.0_real64, 5.0_real64], [n_axes, n_genes])
        gene_to_fam = [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 0]
        family_directions = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT

        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                      family_directions, angular_dispersions, member_counts, status, ierr=ierr)

        call assert_err(ierr, ERR_OK, "family direction exact: ierr")
        call assert_equal_array_int(member_counts, [4, 3, 4], n_families, "family direction exact: member counts")
        call assert_equal_array_int(status, [ERR_OK, STAT_NO_STABLE_DIRECTION, STAT_NO_STABLE_DIRECTION], n_families, &
                                    "family direction exact: status")
        call assert_equal_array_real(family_directions(:, 1), [1.0_real64, 1.0_real64]/sqrt(2.0_real64), n_axes, &
                                     TOL_ULP, "family direction exact: family 1 direction")
        call assert_equal_real(angular_dispersions(1), sqrt(log(2.0_real64)), TOL_ULP, &
                               "family direction exact: family 1 dispersion sqrt(ln 2)")
        call assert_equal_array_real(family_directions(:, 2), [0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                     "family direction exact: family 2 without a stable direction gets the zero vector")
        call assert_equal_real(angular_dispersions(2), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "family direction exact: family 2 dispersion above max is the sentinel")
        call assert_equal_array_real(family_directions(:, 3), [0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                     "family direction exact: family 3 cancelling out gets the zero vector")
        call assert_equal_real(angular_dispersions(3), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "family direction exact: family 3 dispersion is the sentinel")
        ! the assertions above let NaN through; a family whose unit vectors cancel out exactly
        ! must not divide by its zero resultant
        call assert_all_finite(family_directions, n_axes*n_families, "family direction exact: directions")
        call assert_all_finite(angular_dispersions, n_families, "family direction exact: dispersions")
    end subroutine test_family_direction_exact

    !> A zero expression vector is neither counted nor summed: its family's statistics are the
    !| same, bit for bit, as without it.
    subroutine test_family_direction_zero_gene_skipped()
        integer(int32), parameter :: n_axes = 2, n_families = 1
        real(real64) :: with_zero(n_axes, 5), without_zero(n_axes, 4)
        integer(int32) :: gene_to_fam_with(5), gene_to_fam_without(4)
        real(real64) :: directions_with(n_axes, n_families), directions_without(n_axes, n_families)
        real(real64) :: dispersions_with(n_families), dispersions_without(n_families)
        integer(int32) :: counts_with(n_families), counts_without(n_families)
        integer(int32) :: status_with(n_families), status_without(n_families), ierr

        ! family 1 of test_family_direction_exact, once with a zero gene in the middle
        without_zero = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, &
                                0.0_real64, 2.0_real64, 0.0_real64, 0.5_real64], [n_axes, 4])
        with_zero = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
                             0.0_real64, 2.0_real64, 0.0_real64, 0.5_real64], [n_axes, 5])
        gene_to_fam_with = 1
        gene_to_fam_without = 1
        directions_with = UNWRITTEN
        dispersions_with = UNWRITTEN
        counts_with = UNWRITTEN_INT
        status_with = UNWRITTEN_INT

        call compute_family_direction(n_axes, 4, n_families, without_zero, gene_to_fam_without, directions_without, &
                                      dispersions_without, counts_without, status_without, ierr=ierr)
        call assert_err(ierr, ERR_OK, "zero gene: ierr without it")
        call compute_family_direction(n_axes, 5, n_families, with_zero, gene_to_fam_with, directions_with, &
                                      dispersions_with, counts_with, status_with, ierr=ierr)
        call assert_err(ierr, ERR_OK, "zero gene: a zero vector does not abort the call")

        call assert_equal_int(counts_with(1), 4, "zero gene: not counted as a member")
        call assert_equal_int(status_with(1), ERR_OK, "zero gene: family still usable")
        call assert_equal_array_real(directions_with(:, 1), directions_without(:, 1), n_axes, 0.0_real64, &
                                     "zero gene: direction unchanged")
        call assert_equal_real(dispersions_with(1), dispersions_without(1), 0.0_real64, "zero gene: dispersion unchanged")
    end subroutine test_family_direction_zero_gene_skipped

    !> A family needs three members with a direction: none, one, two, or three of which one is
    !| zero, all give STAT_TOO_FEW_MEMBERS, the zero vector and the sentinel.
    subroutine test_family_direction_too_few_members()
        integer(int32), parameter :: n_axes = 2, n_genes = 11, n_families = 5
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families)
        integer(int32) :: member_counts(n_families), status(n_families), ierr, i_family

        ! family 1: one member; family 2: two; family 3: (1,0), (0,1) and a zero vector, so two
        ! with a direction; family 4: none. Family 5 is a usable control: (1,0), (1,1), (0,1).
        expression_vectors = reshape([1.0_real64, 0.0_real64, &
                                      1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                                      1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, &
                                      1.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, &
                                      7.0_real64, 7.0_real64, 8.0_real64, 9.0_real64], [n_axes, n_genes])
        gene_to_fam = [1, 2, 2, 3, 3, 3, 5, 5, 5, 0, 0]
        family_directions = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT

        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                      family_directions, angular_dispersions, member_counts, status, ierr=ierr)

        call assert_err(ierr, ERR_OK, "too few members: ierr")
        call assert_equal_array_int(member_counts, [1, 2, 2, 0, 3], n_families, "too few members: member counts")
        call assert_equal_array_int(status, [STAT_TOO_FEW_MEMBERS, STAT_TOO_FEW_MEMBERS, STAT_TOO_FEW_MEMBERS, &
                                             STAT_TOO_FEW_MEMBERS, ERR_OK], n_families, "too few members: status")
        do i_family = 1, 4
            call assert_equal_array_real(family_directions(:, i_family), [0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                         "too few members: zero vector, never the raw sum, for family "//i_family)
            call assert_equal_real(angular_dispersions(i_family), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                                   "too few members: sentinel dispersion for family "//i_family)
        end do
    end subroutine test_family_direction_too_few_members

    !> Identical and proportional members deviate from their direction by rounding noise only,
    !| within the noise floor: STAT_NO_ANGULAR_VARIATION whatever the minimum, and never NaN.
    subroutine test_family_direction_identical_vectors_never_nan()
        integer(int32), parameter :: n_axes = 3, n_sweep = 200
        real(real64) :: exact(n_axes, 6), sweep(n_axes, 3*n_sweep)
        integer(int32) :: gene_to_fam_exact(6), gene_to_fam_sweep(3*n_sweep)
        real(real64) :: directions_exact(n_axes, 2), dispersions_exact(2)
        real(real64) :: directions_sweep(n_axes, n_sweep), dispersions_sweep(n_sweep)
        integer(int32) :: counts_exact(2), status_exact(2), counts_sweep(n_sweep), status_sweep(n_sweep)
        integer(int32) :: ierr, i_family, i_copy
        real(real64) :: base(n_axes)

        ! family 1: (1,0,0) three times, family 2: (1,0,0), (2,0,0), (5,0,0). Both sum to (3,0,0)
        ! exactly, every deviation is exactly 0, so sigma = 0 -- STAT_NO_ANGULAR_VARIATION even with
        ! a minimum of 0.
        exact = 0.0_real64
        exact(1, :) = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 2.0_real64, 5.0_real64]
        gene_to_fam_exact = [1, 1, 1, 2, 2, 2]
        directions_exact = UNWRITTEN
        dispersions_exact = UNWRITTEN
        call compute_family_direction(n_axes, 6, 2, exact, gene_to_fam_exact, directions_exact, dispersions_exact, &
                                      counts_exact, status_exact, min_angular_dispersion=0.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "identical vectors: ierr")
        call assert_equal_array_int(status_exact, [STAT_NO_ANGULAR_VARIATION, STAT_NO_ANGULAR_VARIATION], 2, &
                                    "identical vectors: sigma 0 is no angular variation, even with min 0")
        call assert_equal_array_real(dispersions_exact, [ANGULAR_DISPERSION_SENTINEL, ANGULAR_DISPERSION_SENTINEL], 2, &
                                     0.0_real64, "identical vectors: sentinel dispersion")
        call assert_equal_array_real(directions_exact(:, 2), [1.0_real64, 0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                     "proportional vectors: direction kept")

        ! 200 families of three copies of a generic vector, scaled by 1, 3 and 1/7: their unit
        ! vectors differ in the last bits, which once gave a rounded R above 1 (sigma NaN), then a
        ! sigma of about 1e-8 with status OK
        do i_family = 1, n_sweep
            base(1) = real(i_family, real64)
            base(2) = 2.0_real64*i_family + 1.0_real64
            base(3) = 3.0_real64*i_family + 7.0_real64
            do i_copy = 1, 3
                gene_to_fam_sweep(3*(i_family - 1) + i_copy) = i_family
            end do
            sweep(:, 3*i_family - 2) = base
            sweep(:, 3*i_family - 1) = 3.0_real64*base
            sweep(:, 3*i_family) = base/7.0_real64
        end do

        ! with the default minimum (0), the noise floor rejects every one of them
        call compute_family_direction(n_axes, 3*n_sweep, n_sweep, sweep, gene_to_fam_sweep, directions_sweep, &
                                      dispersions_sweep, counts_sweep, status_sweep, ierr=ierr)
        call assert_err(ierr, ERR_OK, "identical sweep, default min: ierr")
        call assert_true(all(status_sweep == STAT_NO_ANGULAR_VARIATION), "identical sweep, default min: all no variation")
        call assert_all_finite(dispersions_sweep, n_sweep, "identical sweep, default min: dispersions")

        ! with a minimum of 0 too: the members' deviations are rounding noise, within the noise floor
        call compute_family_direction(n_axes, 3*n_sweep, n_sweep, sweep, gene_to_fam_sweep, directions_sweep, &
                                      dispersions_sweep, counts_sweep, status_sweep, min_angular_dispersion=0.0_real64, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_OK, "identical sweep, min 0: ierr")
        call assert_true(all(status_sweep == STAT_NO_ANGULAR_VARIATION), "identical sweep, min 0: all no variation")
        call assert_equal_array_real(dispersions_sweep, spread(ANGULAR_DISPERSION_SENTINEL, 1, n_sweep), n_sweep, &
                                     0.0_real64, "identical sweep, min 0: all sentinels")
        call assert_all_finite(directions_sweep, n_axes*n_sweep, "identical sweep, min 0: directions")
    end subroutine test_family_direction_identical_vectors_never_nan

    !> The dispersion bounds are inclusive, each reports its own status, and min > max rejects
    !| every dispersion.
    subroutine test_family_direction_dispersion_bounds()
        integer(int32), parameter :: n_axes = 2, n_genes = 7, n_families = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families), sigma_1
        integer(int32) :: member_counts(n_families), status(n_families), ierr

        ! family 1 has sigma = sqrt(ln 2) = 0.8326, family 2 sigma = sqrt(2 ln 3) = 1.4823
        ! (see test_family_direction_exact)
        expression_vectors = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, &
                                      0.0_real64, 2.0_real64, 0.0_real64, 0.5_real64, &
                                      1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64], &
                                     [n_axes, n_genes])
        gene_to_fam = [1, 1, 1, 1, 2, 2, 2]

        ! a maximum of 0.8 rejects family 1 as unstable: no direction
        family_directions = UNWRITTEN
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, max_angular_dispersion=0.8_real64, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_OK, "bounds, max 0.8: ierr")
        call assert_equal_int(status(1), STAT_NO_STABLE_DIRECTION, "bounds, max 0.8: family 1 above max")
        call assert_equal_real(angular_dispersions(1), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "bounds, max 0.8: sentinel")
        call assert_equal_array_real(family_directions(:, 1), [0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                     "bounds, max 0.8: zero vector")

        ! a minimum of 0.9 rejects family 1 as too tight
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, min_angular_dispersion=0.9_real64, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_OK, "bounds, min 0.9: ierr")
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "bounds, min 0.9: family 1 below min")
        call assert_equal_real(angular_dispersions(1), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, "bounds, min 0.9: sentinel")

        ! both bounds exactly at family 1's own sigma: accepted, the bounds are inclusive
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, ierr=ierr)
        sigma_1 = angular_dispersions(1)
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, min_angular_dispersion=sigma_1, &
                                      max_angular_dispersion=sigma_1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "bounds, min = max = sigma: ierr")
        call assert_equal_int(status(1), ERR_OK, "bounds, min = max = sigma: accepted")
        call assert_equal_real(angular_dispersions(1), sigma_1, 0.0_real64, "bounds, min = max = sigma: value")

        ! min 1.0 > max 0.5: family 1 (0.83 < min) is too tight, family 2 (1.48, not below min)
        ! is too loose -- nothing is accepted
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, min_angular_dispersion=1.0_real64, &
                                      max_angular_dispersion=0.5_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "bounds, min > max: ierr")
        call assert_equal_array_int(status, [STAT_NO_ANGULAR_VARIATION, STAT_NO_STABLE_DIRECTION], n_families, &
                                    "bounds, min > max: every family reported")
    end subroutine test_family_direction_dispersion_bounds

    !> Magnitudes whose squares overflow or underflow still give exact directions: the length is
    !| taken on the vector scaled by a power of two that brings its largest component into
    !| [2^-500, 2^500].
    subroutine test_family_direction_extreme_magnitudes()
        integer(int32), parameter :: n_axes = 2, n_genes = 3, n_families = 1
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families), mean_resultant_length
        integer(int32) :: member_counts(n_families), status(n_families), ierr

        ! (1e300, 0), (0, 1e-300) and (1e300, 1e300) -- squares far beyond the range of a double --
        ! have the unit vectors e1, e2 and (1,1)/sqrt(2): the sum is (1 + 1/sqrt(2)) * (1,1),
        ! R = sqrt(2) * (1 + 1/sqrt(2)) / 3 = (sqrt(2) + 1)/3, sigma = sqrt(-2 ln R) = 0.6590, and the
        ! direction is (1,1)/sqrt(2).
        expression_vectors = reshape([1.0e300_real64, 0.0_real64, 0.0_real64, 1.0e-300_real64, &
                                      1.0e300_real64, 1.0e300_real64], [n_axes, n_genes])
        gene_to_fam = 1
        mean_resultant_length = (sqrt(2.0_real64) + 1.0_real64)/3.0_real64
        family_directions = UNWRITTEN

        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                      angular_dispersions, member_counts, status, ierr=ierr)

        call assert_err(ierr, ERR_OK, "extreme magnitudes: ierr")
        call assert_equal_int(status(1), ERR_OK, "extreme magnitudes: status")
        call assert_all_finite(family_directions, n_axes, "extreme magnitudes: direction")
        call assert_equal_array_real(family_directions(:, 1), [1.0_real64, 1.0_real64]/sqrt(2.0_real64), n_axes, &
                                     TOL_ULP, "extreme magnitudes: direction")
        call assert_equal_real(angular_dispersions(1), sqrt(-2.0_real64*log(mean_resultant_length)), TOL_ULP, &
                               "extreme magnitudes: dispersion")
    end subroutine test_family_direction_extreme_magnitudes

    ! ------------------------------------------------------------------------------------------
    ! compute_angular_deviations
    ! ------------------------------------------------------------------------------------------

    !> Angles of exact geometry: parallel 0, orthogonal pi/2, antipodal pi, cos = 1/2 gives pi/3.
    subroutine test_angular_deviations_exact()
        integer(int32), parameter :: n_axes = 2, n_genes = 9, n_families = 3
        real(real64) :: expression_vectors(n_axes, n_genes), family_directions(n_axes, n_families)
        integer(int32) :: gene_to_fam(n_genes), ierr
        real(real64) :: angular_deviations(n_genes)

        ! family 1 points along e1, family 2 has no direction (zero vector), family 3 is e1 at
        ! half length -- the angle does not depend on it
        family_directions = reshape([1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.5_real64, 0.0_real64], &
                                    [n_axes, n_families])
        ! gene 1 (2,0): cos 1, angle 0          gene 2 (0,3): cos 0, angle pi/2
        ! gene 3 (-1,0): cos -1, angle pi       gene 4 (1, sqrt 3): cos 1/2, angle pi/3
        ! gene 5 zero vector: sentinel         gene 6 no family: sentinel
        ! gene 7 in family 2 without direction: sentinel
        ! gene 8 (0,1) against the half-length e1: pi/2
        ! gene 9 (1e300, 1e300): pi/4, although its squared length overflows
        expression_vectors = reshape([2.0_real64, 0.0_real64, 0.0_real64, 3.0_real64, -1.0_real64, 0.0_real64, &
                                      1.0_real64, sqrt(3.0_real64), 0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, &
                                      1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, 1.0e300_real64, 1.0e300_real64], &
                                     [n_axes, n_genes])
        gene_to_fam = [1, 1, 1, 1, 1, 0, 2, 3, 1]
        angular_deviations = UNWRITTEN

        call compute_angular_deviations(n_axes, n_genes, n_families, expression_vectors, family_directions, gene_to_fam, &
                                        angular_deviations, ierr=ierr)

        call assert_err(ierr, ERR_OK, "angular deviations exact: ierr")
        call assert_equal_real(angular_deviations(1), 0.0_real64, 0.0_real64, "angular deviations: parallel is 0")
        call assert_equal_real(angular_deviations(2), PI/2, TOL_ULP, "angular deviations: orthogonal is pi/2")
        call assert_equal_real(angular_deviations(3), PI, TOL_ULP, "angular deviations: antipodal is pi")
        call assert_equal_real(angular_deviations(4), PI/3, TOL_ULP, "angular deviations: cos 1/2 is pi/3")
        call assert_equal_real(angular_deviations(5), ANGULAR_DEVIATIONS_SENTINEL, 0.0_real64, &
                               "angular deviations: zero vector is the sentinel")
        call assert_equal_real(angular_deviations(6), ANGULAR_DEVIATIONS_SENTINEL, 0.0_real64, &
                               "angular deviations: no family is the sentinel")
        call assert_equal_real(angular_deviations(7), ANGULAR_DEVIATIONS_SENTINEL, 0.0_real64, &
                               "angular deviations: family without direction is the sentinel")
        call assert_equal_real(angular_deviations(8), PI/2, TOL_ULP, &
                               "angular deviations: independent of the direction's length")
        call assert_equal_real(angular_deviations(9), PI/4, TOL_ULP, "angular deviations: no overflow for huge values")
    end subroutine test_angular_deviations_exact

    ! ------------------------------------------------------------------------------------------
    ! compute_family_direction_rap
    ! ------------------------------------------------------------------------------------------

    !> Circular means and dispersions, including a family straddling +-pi.
    subroutine test_family_direction_rap_exact()
        integer(int32), parameter :: n_genes = 11, n_families = 3
        real(real64) :: signed_angles(n_genes), family_mean_angles(n_families), angular_dispersions(n_families)
        integer(int32) :: gene_to_fam(n_genes), member_counts(n_families), status(n_families), ierr
        real(real64) :: mean_resultant_length

        ! family 1: 0, pi/2, 0, pi/2 -- unit vectors e1, e2, e1, e2 (cos(pi/2) is 6e-17, not 0):
        !   sum (2,2), R = sqrt(2)/2, sigma = sqrt(ln 2), mean angle pi/4
        ! family 2: pi, -pi + 0.2, pi - 0.2 -- straddling the cut at +-pi. The sines cancel, the
        !   cosines sum to -(1 + 2 cos 0.2): R = (1 + 2 cos 0.2)/3, sigma = sqrt(-2 ln R) = 0.1633,
        !   mean angle pi (or, the same direction, just above -pi)
        ! family 3: 0, pi/2, pi -- the sum is (0, 1) up to the rounding of PI: R = 1/3,
        !   sigma = sqrt(2 ln 3) > max: no stable direction, so no mean angle
        ! gene 11 belongs to no family; its angle must not count
        signed_angles = [0.0_real64, PI/2, 0.0_real64, PI/2, PI, -PI + 0.2_real64, PI - 0.2_real64, &
                         0.0_real64, PI/2, PI, -PI/2]
        gene_to_fam = [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0]
        family_mean_angles = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT

        call compute_family_direction_rap(n_genes, n_families, signed_angles, gene_to_fam, family_mean_angles, &
                                          angular_dispersions, member_counts, status, ierr=ierr)

        call assert_err(ierr, ERR_OK, "rap family exact: ierr")
        call assert_equal_array_int(member_counts, [4, 3, 3], n_families, "rap family exact: member counts")
        call assert_equal_array_int(status, [ERR_OK, ERR_OK, STAT_NO_STABLE_DIRECTION], n_families, &
                                    "rap family exact: status")
        call assert_equal_real(family_mean_angles(1), PI/4, TOL_ULP, "rap family exact: family 1 mean pi/4")
        call assert_equal_real(angular_dispersions(1), sqrt(log(2.0_real64)), TOL_ULP, &
                               "rap family exact: family 1 dispersion sqrt(ln 2)")
        call assert_true(family_mean_angles(2) > -PI .and. family_mean_angles(2) <= PI, &
                         "rap family exact: family 2 mean inside (-pi, pi]")
        call assert_equal_real(abs(wrap_angle(family_mean_angles(2) - PI)), 0.0_real64, TOL_CHAIN, &
                               "rap family exact: family 2 mean is the direction pi")
        mean_resultant_length = (1.0_real64 + 2.0_real64*cos(0.2_real64))/3.0_real64
        call assert_equal_real(angular_dispersions(2), sqrt(-2.0_real64*log(mean_resultant_length)), TOL_CHAIN, &
                               "rap family exact: family 2 dispersion across the cut")
        call assert_equal_real(family_mean_angles(3), FAMILY_MEAN_ANGLE_SENTINEL, 0.0_real64, &
                               "rap family exact: family 3 without a stable direction has the mean sentinel")
        call assert_equal_real(angular_dispersions(3), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "rap family exact: family 3 dispersion above max is the sentinel")
        call assert_all_finite(family_mean_angles, n_families, "rap family exact: means")
        call assert_all_finite(angular_dispersions, n_families, "rap family exact: dispersions")
    end subroutine test_family_direction_rap_exact

    !> One and two members are too few, the mean is the sentinel rather than the raw sum of sines;
    !| identical angles deviate by rounding noise only: no angular variation, never NaN.
    subroutine test_family_direction_rap_too_few_and_clamp()
        integer(int32), parameter :: n_genes = 9, n_families = 4
        real(real64) :: signed_angles(n_genes), family_mean_angles(n_families), angular_dispersions(n_families)
        integer(int32) :: gene_to_fam(n_genes), member_counts(n_families), status(n_families), ierr, i_family

        ! family 1: one member; family 2: two; family 3: 0.3 three times (every deviation from the
        ! mean is rounding noise, so no angular variation even with a minimum of 0); family 4 empty
        signed_angles = [1.0_real64, 1.0_real64, 2.0_real64, 0.3_real64, 0.3_real64, 0.3_real64, &
                         0.5_real64, 0.5_real64, 0.5_real64]
        gene_to_fam = [1, 2, 2, 3, 3, 3, 0, 0, 0]
        family_mean_angles = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT

        call compute_family_direction_rap(n_genes, n_families, signed_angles, gene_to_fam, family_mean_angles, &
                                          angular_dispersions, member_counts, status, min_angular_dispersion=0.0_real64, &
                                          ierr=ierr)

        call assert_err(ierr, ERR_OK, "rap too few: ierr")
        call assert_equal_array_int(member_counts, [1, 2, 3, 0], n_families, "rap too few: member counts")
        call assert_equal_int(status(1), STAT_TOO_FEW_MEMBERS, "rap too few: one member")
        call assert_equal_int(status(2), STAT_TOO_FEW_MEMBERS, "rap too few: two members")
        call assert_equal_int(status(4), STAT_TOO_FEW_MEMBERS, "rap too few: no member")
        do i_family = 1, n_families
            if (i_family == 3) cycle
            call assert_equal_real(family_mean_angles(i_family), FAMILY_MEAN_ANGLE_SENTINEL, 0.0_real64, &
                                   "rap too few: mean angle sentinel, never the raw sum, for family "//i_family)
            call assert_equal_real(angular_dispersions(i_family), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                                   "rap too few: dispersion sentinel for family "//i_family)
        end do
        call assert_all_finite(angular_dispersions, n_families, "rap identical angles")
        call assert_equal_real(family_mean_angles(3), 0.3_real64, TOL_ULP, "rap identical angles: mean kept")
        call assert_equal_int(status(3), STAT_NO_ANGULAR_VARIATION, &
                              "rap identical angles: no angular variation, even with min 0")
        call assert_equal_real(angular_dispersions(3), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "rap identical angles: sentinel dispersion")
    end subroutine test_family_direction_rap_too_few_and_clamp

    ! ------------------------------------------------------------------------------------------
    ! compute_angular_deviations_rap
    ! ------------------------------------------------------------------------------------------

    !> The deviation is the distance around the circle, in [0, pi].
    subroutine test_angular_deviations_rap_wrapping()
        integer(int32), parameter :: n_genes = 7, n_families = 4
        real(real64) :: signed_angles(n_genes), family_mean_angles(n_families), angular_deviations(n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr

        ! family 1 mean pi/2, family 2 mean pi, family 3 mean -pi/2, family 4 has none
        family_mean_angles = [PI/2, PI, -PI/2, FAMILY_MEAN_ANGLE_SENTINEL]
        ! gene 1: -pi/2 against pi/2 is half a turn: |wrap(-pi)| = pi exactly
        ! gene 2: pi/2 against pi/2: 0 exactly
        ! gene 3: -pi + 0.5 against pi: -2 pi + 0.5 wraps to 0.5 -- not 2 pi - 0.5
        ! gene 4: pi against -pi/2: 3 pi/2 wraps to -pi/2, so pi/2
        ! gene 5: pi - 0.25 against pi: 0.25
        ! gene 6: family without a mean: sentinel;  gene 7: no family: sentinel
        signed_angles = [-PI/2, PI/2, -PI + 0.5_real64, PI, PI - 0.25_real64, 1.0_real64, 1.0_real64]
        gene_to_fam = [1, 1, 2, 3, 2, 4, 0]
        angular_deviations = UNWRITTEN

        call compute_angular_deviations_rap(n_genes, n_families, signed_angles, family_mean_angles, gene_to_fam, &
                                            angular_deviations, ierr=ierr)

        call assert_err(ierr, ERR_OK, "rap deviations: ierr")
        call assert_equal_real(angular_deviations(1), PI, 0.0_real64, "rap deviations: half a turn is pi")
        call assert_equal_real(angular_deviations(2), 0.0_real64, 0.0_real64, "rap deviations: same angle is 0")
        call assert_equal_real(angular_deviations(3), 0.5_real64, TOL_ULP, "rap deviations: wraps across -pi")
        call assert_equal_real(angular_deviations(4), PI/2, TOL_ULP, "rap deviations: 3pi/2 wraps to pi/2")
        call assert_equal_real(angular_deviations(5), 0.25_real64, TOL_ULP, "rap deviations: near pi")
        call assert_equal_real(angular_deviations(6), ANGULAR_DEVIATIONS_SENTINEL, 0.0_real64, &
                               "rap deviations: family without mean")
        call assert_equal_real(angular_deviations(7), ANGULAR_DEVIATIONS_SENTINEL, 0.0_real64, "rap deviations: no family")
    end subroutine test_angular_deviations_rap_wrapping

    ! ------------------------------------------------------------------------------------------
    ! compute_relative_angular_deviations
    ! ------------------------------------------------------------------------------------------

    !> Quotients of exactly representable values, and every case that gives the sentinel.
    subroutine test_relative_angular_deviations_exact()
        integer(int32), parameter :: n_genes = 10, n_families = 6
        real(real64) :: angular_deviations(n_genes), angular_dispersions(n_families), relative(n_genes)
        integer(int32) :: gene_to_fam(n_genes), ierr
        real(real64) :: smallest_normal, subnormal

        ! a subnormal dispersion, computed at run time; a flush-to-zero build makes it zero,
        ! which gives the same sentinel
        smallest_normal = tiny(1.0_real64)
        subnormal = smallest_normal/4.0_real64

        ! dispersions: 0.5, 0.25, sentinel, 0, subnormal, smallest normal
        angular_dispersions(1:4) = [0.5_real64, 0.25_real64, ANGULAR_DISPERSION_SENTINEL, 0.0_real64]
        angular_dispersions(5) = subnormal
        angular_dispersions(6) = smallest_normal
        ! gene 1: 1/0.5 = 2          gene 2: 0.5/0.5 = 1        gene 3: 0/0.5 = 0
        ! gene 4: deviation sentinel gene 5: 0.75/0.25 = 3      gene 6: dispersion sentinel
        ! gene 7: no family          gene 8: dispersion 0       gene 9: subnormal dispersion
        ! gene 10: pi / tiny, finite (pi/tiny = 1.4e308 < huge)
        angular_deviations = [1.0_real64, 0.5_real64, 0.0_real64, ANGULAR_DEVIATIONS_SENTINEL, 0.75_real64, &
                              1.0_real64, 1.0_real64, 1.0_real64, PI, PI]
        gene_to_fam = [1, 1, 1, 1, 2, 3, 0, 4, 5, 6]
        relative = UNWRITTEN

        call compute_relative_angular_deviations(n_genes, n_families, angular_deviations, angular_dispersions, &
                                                 gene_to_fam, relative, ierr=ierr)

        call assert_err(ierr, ERR_OK, "relative: ierr")
        call assert_equal_array_real(relative(1:9), [2.0_real64, 1.0_real64, 0.0_real64, RELATIVE_SENTINEL, 3.0_real64, &
                                     RELATIVE_SENTINEL, RELATIVE_SENTINEL, RELATIVE_SENTINEL, RELATIVE_SENTINEL], 9, &
                                     0.0_real64, "relative: exact quotients and sentinels")
        call assert_true(ieee_is_finite(relative(10)), "relative: pi over the smallest normal dispersion is finite")
        call assert_equal_real(relative(10)/(PI/smallest_normal), 1.0_real64, TOL_ULP, "relative: pi / tiny")
    end subroutine test_relative_angular_deviations_exact

    ! ------------------------------------------------------------------------------------------
    ! compute_angle_outlier_threshold and flag_angle_outliers
    ! ------------------------------------------------------------------------------------------

    !> The type-7 quantile of the values that exist, in both tiers.
    subroutine test_threshold_quantile()
        integer(int32), parameter :: n_genes = 5
        real(real64) :: values(n_genes), threshold
        integer(int32) :: perm(n_genes), ierr

        ! valid values sorted: 0.5, 1, 2, 4 (n = 4); the sentinel does not count.
        ! Type 7: rank h = level*(n - 1) + 1, value = x(floor h) + frac(h)*(x(floor h + 1) - x(floor h))
        values = [0.5_real64, 4.0_real64, RELATIVE_SENTINEL, 1.0_real64, 2.0_real64]

        ! level 0.5: h = 2.5, so 1 + 0.5*(2 - 1) = 1.5
        threshold = UNWRITTEN
        call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=0.5_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "threshold 0.5: ierr")
        call assert_equal_real(threshold, 1.5_real64, 0.0_real64, "threshold 0.5: 1.5")
        ! level 0: h = 1, the smallest value; level 1: h = 4, the largest
        call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=0.0_real64, ierr=ierr)
        call assert_equal_real(threshold, 0.5_real64, 0.0_real64, "threshold 0: smallest value")
        call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=1.0_real64, ierr=ierr)
        call assert_equal_real(threshold, 4.0_real64, 0.0_real64, "threshold 1: largest value")
        ! the default 0.95: h = 3.85, so 2 + 0.85*(4 - 2) = 3.7 (0.95 is not exact in binary)
        call compute_angle_outlier_threshold(n_genes, values, threshold, ierr=ierr)
        call assert_err(ierr, ERR_OK, "threshold default: ierr")
        call assert_equal_real(threshold, 3.7_real64, TOL_ULP*4, "threshold default 0.95: 3.7")

        ! the expert tier takes the work array and leaves the sorted valid indices in front
        threshold = UNWRITTEN
        perm = UNWRITTEN_INT
        call compute_angle_outlier_threshold_expert(n_genes, values, perm, threshold, quantile_level=0.5_real64, &
                                                    ierr=ierr)
        call assert_err(ierr, ERR_OK, "threshold expert: ierr")
        call assert_equal_real(threshold, 1.5_real64, 0.0_real64, "threshold expert: 1.5")
        call assert_equal_array_int(perm, [1, 4, 5, 2, 0], n_genes, "threshold expert: sorted valid indices, then 0")
    end subroutine test_threshold_quantile

    !> Fewer than two values leave nothing to rank: the threshold is huge and nothing is flagged.
    subroutine test_threshold_fewer_than_two_values()
        real(real64) :: one_valid(3), none_valid(2), single(1), threshold
        logical(c_bool) :: is_outlier(3)
        integer(int32) :: ierr

        one_valid = [RELATIVE_SENTINEL, 3.0_real64, RELATIVE_SENTINEL]
        none_valid = RELATIVE_SENTINEL
        single = 0.7_real64

        threshold = UNWRITTEN
        call compute_angle_outlier_threshold(3, one_valid, threshold, ierr=ierr)
        call assert_err(ierr, ERR_OK, "one valid value: ierr")
        call assert_equal_real(threshold, huge(1.0_real64), 0.0_real64, "one valid value: huge")
        is_outlier = .true.
        call flag_angle_outliers(3, one_valid, threshold, is_outlier, ierr=ierr)
        call assert_false(any(is_outlier), "one valid value: nothing flagged")

        threshold = UNWRITTEN
        call compute_angle_outlier_threshold(2, none_valid, threshold, quantile_level=0.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "no valid value: ierr")
        call assert_equal_real(threshold, huge(1.0_real64), 0.0_real64, "no valid value: huge")

        threshold = UNWRITTEN
        call compute_angle_outlier_threshold(1, single, threshold, ierr=ierr)
        call assert_err(ierr, ERR_OK, "single gene: ierr")
        call assert_equal_real(threshold, huge(1.0_real64), 0.0_real64, "single gene: huge")
    end subroutine test_threshold_fewer_than_two_values

    !> Ties all reach the threshold; zeros never count as outliers, even when they are the
    !| threshold.
    subroutine test_threshold_and_flags_with_ties_and_zeros()
        real(real64) :: ties(4), zeros(3), threshold
        logical(c_bool) :: flags_ties(4), flags_zeros(3)
        integer(int32) :: ierr

        ! 2, 2, 2: every quantile is 2, and no value lies strictly above it -- all equal, none flagged
        ties = [2.0_real64, RELATIVE_SENTINEL, 2.0_real64, 2.0_real64]
        call compute_angle_outlier_threshold(4, ties, threshold, quantile_level=0.5_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "ties: ierr")
        call assert_equal_real(threshold, 2.0_real64, 0.0_real64, "ties: threshold 2")
        flags_ties = .true.
        call flag_angle_outliers(4, ties, threshold, flags_ties, ierr=ierr)
        call assert_false(any(flags_ties), "ties: values equal to the threshold are not flagged")

        ! 0, 0: the threshold is 0, but a deviation of zero is never an outlier
        zeros = [0.0_real64, 0.0_real64, RELATIVE_SENTINEL]
        call compute_angle_outlier_threshold(3, zeros, threshold, ierr=ierr)
        call assert_equal_real(threshold, 0.0_real64, 0.0_real64, "zeros: threshold 0")
        flags_zeros = .true.
        call flag_angle_outliers(3, zeros, threshold, flags_zeros, ierr=ierr)
        call assert_false(any(flags_zeros), "zeros: nothing flagged")
    end subroutine test_threshold_and_flags_with_ties_and_zeros

    !> is_outlier = v > 0 .and. v > threshold, whatever threshold the caller passes.
    subroutine test_flag_rule()
        real(real64) :: values(4)
        logical(c_bool) :: is_outlier(4)
        integer(int32) :: ierr

        values = [RELATIVE_SENTINEL, 0.0_real64, 0.5_real64, 3.0_real64]

        ! a negative threshold must not reach the sentinel or the zero
        is_outlier = .true.
        call flag_angle_outliers(4, values, -5.0_real64, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_OK, "flag rule: ierr")
        call assert_equal_array_logical(is_outlier, [.false._c_bool, .false._c_bool, .true._c_bool, .true._c_bool], 4, &
                                        "flag rule: negative threshold")
        ! a value exactly equal to the threshold is not flagged
        call flag_angle_outliers(4, values, 0.5_real64, is_outlier, ierr=ierr)
        call assert_equal_array_logical(is_outlier, [.false._c_bool, .false._c_bool, .false._c_bool, .true._c_bool], 4, &
                                        "flag rule: a value equal to the threshold is not flagged")
        call flag_angle_outliers(4, values, 0.6_real64, is_outlier, ierr=ierr)
        call assert_equal_array_logical(is_outlier, [.false._c_bool, .false._c_bool, .false._c_bool, .true._c_bool], 4, &
                                        "flag rule: above one value")
        call flag_angle_outliers(4, values, huge(1.0_real64), is_outlier, ierr=ierr)
        call assert_false(any(is_outlier), "flag rule: huge threshold flags nothing")
    end subroutine test_flag_rule

    !> A ranking cut with a strict comparison: the largest value is flagged exactly when it lies
    !| strictly above the threshold -- at every level for a unique maximum below level 1, never for
    !| top values that tie at the threshold.
    subroutine test_largest_flagged_iff_above_threshold()
        integer(int32), parameter :: n_genes = 24
        real(real64) :: values(n_genes), threshold, level, tied_top(5)
        logical(c_bool) :: is_outlier(n_genes), tied_flags(5)
        integer(int32) :: ierr, i_gene, i_level, i_largest

        ! distinct values 0.5 * (7 i mod 23), with every fifth gene a sentinel
        do i_gene = 1, n_genes
            values(i_gene) = 0.5_real64*real(mod(7*i_gene, 23), real64)
            if (mod(i_gene, 5) == 0) values(i_gene) = RELATIVE_SENTINEL
        end do
        i_largest = maxloc(values, dim=1)

        do i_level = 0, 10
            level = real(i_level, real64)/10.0_real64
            call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=level, ierr=ierr)
            call assert_err(ierr, ERR_OK, "largest: ierr at level "//level)
            call flag_angle_outliers(n_genes, values, threshold, is_outlier, ierr=ierr)
            call assert_true(is_outlier(i_largest) .eqv. (values(i_largest) > threshold), &
                             "largest: flagged iff strictly above the threshold, at level "//level)
        end do

        ! a unique maximum above the interpolated threshold: level 0.9 over the 20 values puts the
        ! type-7 rank at 0.9 * 19 + 1 = 18.1, below the two largest, so the maximum is flagged
        call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=0.9_real64, ierr=ierr)
        call flag_angle_outliers(n_genes, values, threshold, is_outlier, ierr=ierr)
        call assert_true(is_outlier(i_largest), "largest: a unique maximum above the threshold is flagged")
        ! level 1 makes the threshold the maximum itself, which is then not above it
        call compute_angle_outlier_threshold(n_genes, values, threshold, quantile_level=1.0_real64, ierr=ierr)
        call flag_angle_outliers(n_genes, values, threshold, is_outlier, ierr=ierr)
        call assert_false(any(is_outlier), "largest: at level 1 nothing lies above the threshold")

        ! 1, 2, 5, 5, 5: level 0.5 has rank 3, the first of the tied top values, so the threshold is
        ! 5 and none of the tied maxima is flagged
        tied_top = [1.0_real64, 5.0_real64, 2.0_real64, 5.0_real64, 5.0_real64]
        call compute_angle_outlier_threshold(5, tied_top, threshold, quantile_level=0.5_real64, ierr=ierr)
        call assert_equal_real(threshold, 5.0_real64, 0.0_real64, "tied top: threshold equals the tied maxima")
        tied_flags = .true.
        call flag_angle_outliers(5, tied_top, threshold, tied_flags, ierr=ierr)
        call assert_false(any(tied_flags), "tied top: maxima tied at the threshold are not flagged")
    end subroutine test_largest_flagged_iff_above_threshold

    ! ------------------------------------------------------------------------------------------
    ! detect_angle_outliers
    ! ------------------------------------------------------------------------------------------

    !> The spherical pipeline end to end: every family status and gene status, hand-derived
    !| relative angular deviations, the flags, the expert tier, and every output of the steps
    !| called one by one.
    subroutine test_detect_angle_outliers()
        integer(int32), parameter :: n_axes = 2, n_genes = 17, n_families = 6
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families)
        real(real64) :: relative(n_genes), threshold
        integer(int32) :: member_counts(n_families), status(n_families), gene_status(n_genes), ierr
        logical(c_bool) :: is_outlier(n_genes), expected_outlier(n_genes)
        ! the expert tier and the steps, for comparison
        real(real64) :: expert_directions(n_axes, n_families), expert_dispersions(n_families), expert_relative(n_genes)
        real(real64) :: expert_threshold, tmp_angular_deviations(n_genes)
        integer(int32) :: expert_counts(n_families), expert_status(n_families), expert_gene_status(n_genes)
        integer(int32) :: tmp_perm(n_genes)
        logical(c_bool) :: expert_is_outlier(n_genes)
        real(real64) :: step_directions(n_axes, n_families), step_dispersions(n_families), step_deviations(n_genes)
        real(real64) :: step_relative(n_genes)
        integer(int32) :: step_counts(n_families), step_status(n_families)
        real(real64) :: level, relative_family_1, relative_family_5, sigma_family_5, expected_relative(n_genes)

        ! family 1, genes 1-4: e1, e1, e2, e2 (scaled) -- direction (1,1)/sqrt(2), sigma sqrt(ln 2),
        !   every angle pi/4, so each relative angular deviation is (pi/4)/sqrt(ln 2) = 0.9433
        ! gene 5: zero vector in family 1           -> STAT_ZERO_VECTOR, not counted
        ! family 2, genes 6-7: two members          -> STAT_TOO_FEW_MEMBERS
        ! family 3, genes 8-10: e1, e2, -e1, R = 1/3 -> STAT_NO_STABLE_DIRECTION
        ! gene 11: no family                        -> STAT_NO_FAMILY
        ! family 5, genes 12-14: e1, (1,1), e2 -- R = (sqrt(2) + 1)/3, sigma = 0.6590, direction
        !   (1,1)/sqrt(2); angles pi/4, 0, pi/4, relative (pi/4)/sigma = 1.1918, 0, 1.1918
        ! family 4, genes 15-17: e1 three times     -> sigma 0, STAT_NO_ANGULAR_VARIATION
        ! family 6: no gene at all                  -> STAT_TOO_FEW_MEMBERS, reported although
        !   it is the last family
        expression_vectors = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, &
                                      0.0_real64, 0.5_real64, 0.0_real64, 0.0_real64, &
                                      1.0_real64, 2.0_real64, 2.0_real64, 1.0_real64, &
                                      1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64, &
                                      4.0_real64, 4.0_real64, &
                                      1.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, &
                                      1.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 1.0_real64, 0.0_real64], &
                                     [n_axes, n_genes])
        gene_to_fam = [1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 0, 5, 5, 5, 4, 4, 4]
        relative_family_1 = (PI/4)/sqrt(log(2.0_real64))
        sigma_family_5 = sqrt(-2.0_real64*log((sqrt(2.0_real64) + 1.0_real64)/3.0_real64))
        relative_family_5 = (PI/4)/sigma_family_5

        ! 7 relative angular deviations exist: 0, 0.9433 (4x), 1.1918 (2x). Level 0.75 puts the
        ! type-7 rank at 0.75 * 6 + 1 = 5.5 (exact), halfway between 0.9433 and 1.1918 -- so exactly
        ! the two 1.1918 genes are flagged, however their last bits differ
        level = 0.75_real64
        expected_outlier = .false.
        expected_outlier([12, 14]) = .true.

        family_directions = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT
        relative = UNWRITTEN
        threshold = UNWRITTEN
        is_outlier = .true.
        gene_status = UNWRITTEN_INT

        call detect_angle_outliers(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                   angular_dispersions, member_counts, status, relative, threshold, is_outlier, &
                                   gene_status, quantile_level=level, ierr=ierr)

        call assert_err(ierr, ERR_OK, "detect: ierr")
        call assert_equal_array_int(member_counts, [4, 2, 3, 3, 3, 0], n_families, "detect: member counts")
        call assert_equal_array_int(status, [ERR_OK, STAT_TOO_FEW_MEMBERS, STAT_NO_STABLE_DIRECTION, &
                                             STAT_NO_ANGULAR_VARIATION, ERR_OK, STAT_TOO_FEW_MEMBERS], n_families, &
                                    "detect: family status")
        call assert_equal_array_int(gene_status, [ERR_OK, ERR_OK, ERR_OK, ERR_OK, STAT_ZERO_VECTOR, &
                                                  STAT_TOO_FEW_MEMBERS, STAT_TOO_FEW_MEMBERS, &
                                                  STAT_NO_STABLE_DIRECTION, STAT_NO_STABLE_DIRECTION, &
                                                  STAT_NO_STABLE_DIRECTION, STAT_NO_FAMILY, ERR_OK, ERR_OK, ERR_OK, &
                                                  STAT_NO_ANGULAR_VARIATION, STAT_NO_ANGULAR_VARIATION, &
                                                  STAT_NO_ANGULAR_VARIATION], n_genes, "detect: gene status")
        ! the sentinel exactly where the gene status is not zero
        expected_relative = RELATIVE_SENTINEL
        expected_relative(1:4) = relative_family_1
        expected_relative(12) = relative_family_5
        expected_relative(13) = 0.0_real64
        expected_relative(14) = relative_family_5
        call assert_equal_array_real(relative, expected_relative, n_genes, TOL_ULP*4, "detect: relative deviations")
        call assert_equal_real(angular_dispersions(5), sigma_family_5, TOL_ULP, "detect: family 5 dispersion")
        call assert_equal_array_real(family_directions(:, 5), [1.0_real64, 1.0_real64]/sqrt(2.0_real64), n_axes, &
                                     TOL_ULP, "detect: family 5 direction")
        call assert_equal_array_real(family_directions(:, 6), [0.0_real64, 0.0_real64], n_axes, 0.0_real64, &
                                     "detect: empty family has the zero vector")
        call assert_in_range_real(threshold, relative_family_1, relative_family_5 + TOL_ULP*4, &
                                  "detect: threshold between the two groups")
        call assert_equal_array_logical(is_outlier, expected_outlier, n_genes, "detect: exactly the two largest flagged")

        ! the expert tier with its own work arrays gives the same, bit for bit
        tmp_angular_deviations = UNWRITTEN
        tmp_perm = UNWRITTEN_INT
        call detect_angle_outliers_expert(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, &
                                          expert_directions, expert_dispersions, expert_counts, expert_status, &
                                          expert_relative, expert_threshold, expert_is_outlier, expert_gene_status, &
                                          tmp_angular_deviations, tmp_perm, quantile_level=level, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect expert: ierr")
        call assert_equal_array_real(expert_relative, relative, n_genes, 0.0_real64, "detect expert: relative")
        call assert_equal_real(expert_threshold, threshold, 0.0_real64, "detect expert: threshold")
        call assert_equal_array_logical(expert_is_outlier, is_outlier, n_genes, "detect expert: flags")
        call assert_equal_array_int(expert_gene_status, gene_status, n_genes, "detect expert: gene status")
        call assert_equal_array_int(expert_status, status, n_families, "detect expert: status")
        call assert_equal_array_int(expert_counts, member_counts, n_families, "detect expert: counts")
        call assert_equal_array_real(expert_directions, family_directions, n_axes*n_families, 0.0_real64, &
                                     "detect expert: directions")
        call assert_equal_array_real(expert_dispersions, angular_dispersions, n_families, 0.0_real64, &
                                     "detect expert: dispersions")
        call assert_equal_array_real(tmp_angular_deviations(1:4), [PI/4, PI/4, PI/4, PI/4], 4, TOL_ULP, &
                                     "detect expert: the work array held the angles")

        ! and the same as the steps, called one by one
        call compute_family_direction(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, step_directions, &
                                      step_dispersions, step_counts, step_status, ierr=ierr)
        call compute_angular_deviations(n_axes, n_genes, n_families, expression_vectors, step_directions, gene_to_fam, &
                                        step_deviations, ierr=ierr)
        call compute_relative_angular_deviations(n_genes, n_families, step_deviations, step_dispersions, gene_to_fam, &
                                                 step_relative, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect steps: ierr")
        call assert_equal_array_real(step_directions, family_directions, n_axes*n_families, 0.0_real64, &
                                     "detect: directions as the steps give them")
        call assert_equal_array_real(step_dispersions, angular_dispersions, n_families, 0.0_real64, &
                                     "detect: dispersions as the steps give them")
        call assert_equal_array_int(step_counts, member_counts, n_families, "detect: counts as the steps give them")
        call assert_equal_array_int(step_status, status, n_families, "detect: status as the steps give it")
        call assert_equal_array_real(step_deviations, tmp_angular_deviations, n_genes, 0.0_real64, &
                                     "detect: angles as the steps give them")
        call assert_equal_array_real(step_relative, relative, n_genes, 0.0_real64, &
                                     "detect: relative deviations as the steps give them")

        ! the default level is 0.95
        call detect_angle_outliers(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                   angular_dispersions, member_counts, status, relative, threshold, is_outlier, &
                                   gene_status, ierr=ierr)
        call compute_angle_outlier_threshold(n_genes, relative, expert_threshold, quantile_level=0.95_real64, ierr=ierr)
        call assert_equal_real(threshold, expert_threshold, 0.0_real64, "detect: the default level is 0.95")
    end subroutine test_detect_angle_outliers

    !> No family has statistics: the threshold is huge and nothing is flagged, and nothing aborts.
    subroutine test_detect_angle_outliers_no_family_has_statistics()
        integer(int32), parameter :: n_axes = 2, n_genes = 4, n_families = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        integer(int32) :: gene_to_fam(n_genes)
        real(real64) :: family_directions(n_axes, n_families), angular_dispersions(n_families)
        real(real64) :: relative(n_genes), threshold
        integer(int32) :: member_counts(n_families), status(n_families), gene_status(n_genes), ierr
        logical(c_bool) :: is_outlier(n_genes)

        ! two families of two members each
        expression_vectors = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                                      1.0_real64, 1.0_real64, 2.0_real64, 1.0_real64], [n_axes, n_genes])
        gene_to_fam = [1, 1, 2, 2]
        relative = UNWRITTEN
        threshold = UNWRITTEN
        is_outlier = .true.

        call detect_angle_outliers(n_axes, n_genes, n_families, expression_vectors, gene_to_fam, family_directions, &
                                   angular_dispersions, member_counts, status, relative, threshold, is_outlier, &
                                   gene_status, ierr=ierr)

        call assert_err(ierr, ERR_OK, "no statistics: ierr")
        call assert_true(all(status == STAT_TOO_FEW_MEMBERS), "no statistics: every family too small")
        call assert_true(all(gene_status == STAT_TOO_FEW_MEMBERS), "no statistics: every gene reports its family")
        call assert_equal_array_real(relative, spread(RELATIVE_SENTINEL, 1, n_genes), n_genes, 0.0_real64, &
                                     "no statistics: all relative deviations are sentinels")
        call assert_equal_real(threshold, huge(1.0_real64), 0.0_real64, "no statistics: threshold huge")
        call assert_false(any(is_outlier), "no statistics: nothing flagged")
    end subroutine test_detect_angle_outliers_no_family_has_statistics

    ! ------------------------------------------------------------------------------------------
    ! detect_angle_outliers_rap
    ! ------------------------------------------------------------------------------------------

    !> The RAP pipeline end to end, across the +-pi cut, with the expert tier and every output of the
    !| steps called one by one.
    subroutine test_detect_angle_outliers_rap()
        integer(int32), parameter :: n_genes = 10, n_families = 3
        real(real64) :: signed_angles(n_genes), family_mean_angles(n_families), angular_dispersions(n_families)
        real(real64) :: relative(n_genes), threshold
        integer(int32) :: gene_to_fam(n_genes), member_counts(n_families), status(n_families), gene_status(n_genes)
        integer(int32) :: ierr
        logical(c_bool) :: is_outlier(n_genes), expected_outlier(n_genes)
        real(real64) :: expert_means(n_families), expert_dispersions(n_families), expert_relative(n_genes)
        real(real64) :: expert_threshold, tmp_angular_deviations(n_genes)
        integer(int32) :: expert_counts(n_families), expert_status(n_families), expert_gene_status(n_genes)
        integer(int32) :: tmp_perm(n_genes)
        logical(c_bool) :: expert_is_outlier(n_genes)
        real(real64) :: step_means(n_families), step_dispersions(n_families), step_deviations(n_genes)
        real(real64) :: step_relative(n_genes)
        integer(int32) :: step_counts(n_families), step_status(n_families)
        real(real64) :: level, sigma_family_2, relative_family_1, relative_family_2, expected_relative(n_genes)

        ! family 1, genes 1-4: 0, pi/2, 0, pi/2 -- mean pi/4, sigma sqrt(ln 2), each deviation
        !   pi/4, relative (pi/4)/sqrt(ln 2) = 0.9433
        ! family 2, genes 5-7: pi, -pi + 0.2, pi - 0.2 -- mean pi, sigma = 0.1633 (see
        !   test_family_direction_rap_exact); deviations 0, 0.2, 0.2 (around the circle, not
        !   2 pi - 0.2), relative 0, 1.2247, 1.2247
        ! family 3, genes 8-9: too few; gene 10: no family
        signed_angles = [0.0_real64, PI/2, 0.0_real64, PI/2, PI, -PI + 0.2_real64, PI - 0.2_real64, &
                         1.0_real64, 2.0_real64, 3.0_real64]
        gene_to_fam = [1, 1, 1, 1, 2, 2, 2, 3, 3, 0]
        sigma_family_2 = sqrt(-2.0_real64*log((1.0_real64 + 2.0_real64*cos(0.2_real64))/3.0_real64))
        relative_family_1 = (PI/4)/sqrt(log(2.0_real64))
        relative_family_2 = 0.2_real64/sigma_family_2

        ! 7 values: 0, 0.9433 (4x), 1.2247 (2x): level 0.75 puts the threshold halfway between
        ! 0.9433 and 1.2247, so exactly the two largest are flagged (see test_detect_angle_outliers)
        level = 0.75_real64
        expected_outlier = .false.
        expected_outlier([6, 7]) = .true.

        family_mean_angles = UNWRITTEN
        angular_dispersions = UNWRITTEN
        member_counts = UNWRITTEN_INT
        status = UNWRITTEN_INT
        relative = UNWRITTEN
        threshold = UNWRITTEN
        is_outlier = .true.
        gene_status = UNWRITTEN_INT

        call detect_angle_outliers_rap(n_genes, n_families, signed_angles, gene_to_fam, family_mean_angles, &
                                       angular_dispersions, member_counts, status, relative, threshold, is_outlier, &
                                       gene_status, quantile_level=level, ierr=ierr)

        call assert_err(ierr, ERR_OK, "detect rap: ierr")
        call assert_equal_array_int(member_counts, [4, 3, 2], n_families, "detect rap: member counts")
        call assert_equal_array_int(status, [ERR_OK, ERR_OK, STAT_TOO_FEW_MEMBERS], n_families, "detect rap: status")
        call assert_equal_array_int(gene_status, [ERR_OK, ERR_OK, ERR_OK, ERR_OK, ERR_OK, ERR_OK, ERR_OK, &
                                                  STAT_TOO_FEW_MEMBERS, STAT_TOO_FEW_MEMBERS, STAT_NO_FAMILY], n_genes, &
                                    "detect rap: gene status")
        call assert_equal_real(family_mean_angles(3), FAMILY_MEAN_ANGLE_SENTINEL, 0.0_real64, &
                               "detect rap: too small family has the mean sentinel")
        ! family 1 at (pi/4)/sqrt(ln 2); the gene at pi deviates by 0 and the two across the cut by
        ! 0.2; the sentinel where the gene status is not 0
        expected_relative = RELATIVE_SENTINEL
        expected_relative(1:4) = relative_family_1
        expected_relative(5) = 0.0_real64
        expected_relative(6:7) = relative_family_2
        call assert_equal_array_real(relative, expected_relative, n_genes, TOL_CHAIN, "detect rap: relative deviations")
        call assert_equal_array_logical(is_outlier, expected_outlier, n_genes, "detect rap: exactly the two largest")

        tmp_angular_deviations = UNWRITTEN
        tmp_perm = UNWRITTEN_INT
        call detect_angle_outliers_rap_expert(n_genes, n_families, signed_angles, gene_to_fam, expert_means, &
                                              expert_dispersions, expert_counts, expert_status, expert_relative, &
                                              expert_threshold, expert_is_outlier, expert_gene_status, &
                                              tmp_angular_deviations, tmp_perm, quantile_level=level, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect rap expert: ierr")
        call assert_equal_array_real(expert_relative, relative, n_genes, 0.0_real64, "detect rap expert: relative")
        call assert_equal_real(expert_threshold, threshold, 0.0_real64, "detect rap expert: threshold")
        call assert_equal_array_logical(expert_is_outlier, is_outlier, n_genes, "detect rap expert: flags")
        call assert_equal_array_int(expert_gene_status, gene_status, n_genes, "detect rap expert: gene status")
        call assert_equal_array_int(expert_status, status, n_families, "detect rap expert: status")
        call assert_equal_array_int(expert_counts, member_counts, n_families, "detect rap expert: counts")
        call assert_equal_array_real(expert_means, family_mean_angles, n_families, 0.0_real64, "detect rap expert: means")
        call assert_equal_array_real(expert_dispersions, angular_dispersions, n_families, 0.0_real64, &
                                     "detect rap expert: dispersions")

        call compute_family_direction_rap(n_genes, n_families, signed_angles, gene_to_fam, step_means, step_dispersions, &
                                          step_counts, step_status, ierr=ierr)
        call compute_angular_deviations_rap(n_genes, n_families, signed_angles, step_means, gene_to_fam, &
                                            step_deviations, ierr=ierr)
        call compute_relative_angular_deviations(n_genes, n_families, step_deviations, step_dispersions, gene_to_fam, &
                                                 step_relative, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect rap steps: ierr")
        call assert_equal_array_real(step_means, family_mean_angles, n_families, 0.0_real64, &
                                     "detect rap: means as the steps give them")
        call assert_equal_array_real(step_dispersions, angular_dispersions, n_families, 0.0_real64, &
                                     "detect rap: dispersions as the steps give them")
        call assert_equal_array_int(step_counts, member_counts, n_families, "detect rap: counts as the steps give them")
        call assert_equal_array_int(step_status, status, n_families, "detect rap: status as the steps give it")
        call assert_equal_array_real(step_deviations, tmp_angular_deviations, n_genes, 0.0_real64, &
                                     "detect rap: deviations as the steps give them")
        call assert_equal_array_real(step_relative, relative, n_genes, 0.0_real64, &
                                     "detect rap: relative deviations as the steps give them")
    end subroutine test_detect_angle_outliers_rap

    ! ------------------------------------------------------------------------------------------
    ! validation: every check on its own, with its code and the argument it blames
    ! ------------------------------------------------------------------------------------------

    subroutine test_validation_compute_family_direction()
        real(real64) :: expression_vectors(2, 3), bad(2, 3), directions(2, 1), dispersions(1)
        integer(int32) :: gene_to_fam(3), counts(1), status(1), ierr

        expression_vectors = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64], [2, 3])
        gene_to_fam = [1, 1, 0]

        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_OK, "family direction: valid input, the sentinel 0 included")
        call compute_family_direction(0, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "family direction: n_axes 0", 1)
        call compute_family_direction(-1, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: n_axes negative", 1)
        call compute_family_direction(2, 0, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "family direction: n_genes 0", 2)
        call compute_family_direction(2, 3, 0, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "family direction: n_families 0", 3)
        bad = expression_vectors
        bad(2, 2) = NAN_VALUE
        call compute_family_direction(2, 3, 1, bad, gene_to_fam, directions, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "family direction: NaN expression", 4)
        bad(2, 2) = -INF_VALUE
        call compute_family_direction(2, 3, 1, bad, gene_to_fam, directions, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "family direction: infinite expression", 4)
        call compute_family_direction(2, 3, 1, expression_vectors, [1, -1, 1], directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: gene_to_fam below 1", 5)
        call compute_family_direction(2, 3, 1, expression_vectors, [1, 2, 1], directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: gene_to_fam above n_families", 5)
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      min_angular_dispersion=-0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: negative min", 10)
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      min_angular_dispersion=NAN_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "family direction: NaN min", 10)
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      max_angular_dispersion=-0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: negative max", 11)
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      max_angular_dispersion=INF_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "family direction: infinite max", 11)
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      max_angular_dispersion=MAX_DISPERSION_LIMIT, ierr=ierr)
        call assert_err(ierr, ERR_OK, "family direction: max at the limit 5")
        call compute_family_direction(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                      max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "family direction: max above the limit 5", 11)
    end subroutine test_validation_compute_family_direction

    subroutine test_validation_compute_angular_deviations()
        real(real64) :: expression_vectors(2, 2), directions(2, 1), bad_directions(2, 1), bad(2, 2), deviations(2)
        integer(int32) :: gene_to_fam(2), ierr

        expression_vectors = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64], [2, 2])
        directions = reshape([-1.0_real64, 1.0_real64], [2, 1])
        gene_to_fam = [1, 0]

        call compute_angular_deviations(2, 2, 1, expression_vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_OK, "angular deviations: valid input, bounds -1 and 1 included")
        call compute_angular_deviations(0, 2, 1, expression_vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "angular deviations: n_axes 0", 1)
        call compute_angular_deviations(2, 0, 1, expression_vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "angular deviations: n_genes 0", 2)
        call compute_angular_deviations(2, 2, 0, expression_vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "angular deviations: n_families 0", 3)
        bad = expression_vectors
        bad(1, 1) = NAN_VALUE
        call compute_angular_deviations(2, 2, 1, bad, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "angular deviations: NaN expression", 4)
        bad(1, 1) = INF_VALUE
        call compute_angular_deviations(2, 2, 1, bad, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "angular deviations: infinite expression", 4)
        bad_directions = reshape([1.5_real64, 0.0_real64], [2, 1])
        call compute_angular_deviations(2, 2, 1, expression_vectors, bad_directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "angular deviations: direction above 1", 5)
        bad_directions = reshape([0.0_real64, -1.5_real64], [2, 1])
        call compute_angular_deviations(2, 2, 1, expression_vectors, bad_directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "angular deviations: direction below -1", 5)
        bad_directions = reshape([0.0_real64, NAN_VALUE], [2, 1])
        call compute_angular_deviations(2, 2, 1, expression_vectors, bad_directions, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "angular deviations: NaN direction", 5)
        call compute_angular_deviations(2, 2, 1, expression_vectors, directions, [2, 1], deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "angular deviations: gene_to_fam above n_families", 6)
        call compute_angular_deviations(2, 2, 1, expression_vectors, directions, [-3, 1], deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "angular deviations: gene_to_fam below 1", 6)
    end subroutine test_validation_compute_angular_deviations

    subroutine test_validation_compute_family_direction_rap()
        real(real64) :: signed_angles(3), means(1), dispersions(1)
        integer(int32) :: gene_to_fam(3), counts(1), status(1), ierr

        ! pi itself is accepted, -pi is not: the range is (-pi, pi]
        signed_angles = [PI, 0.0_real64, nearest(-PI, 1.0_real64)]
        gene_to_fam = [1, 1, 0]

        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "rap family: valid input, pi and just above -pi included")
        call compute_family_direction_rap(0, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "rap family: n_genes 0", 1)
        call compute_family_direction_rap(3, 0, signed_angles, gene_to_fam, means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "rap family: n_families 0", 2)
        call compute_family_direction_rap(3, -2, signed_angles, gene_to_fam, means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: n_families negative", 2)
        call compute_family_direction_rap(3, 1, [0.0_real64, -PI, 0.0_real64], gene_to_fam, means, dispersions, counts, &
                                          status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: -pi is outside (-pi, pi]", 3)
        call compute_family_direction_rap(3, 1, [0.0_real64, nearest(PI, 1.0_real64), 0.0_real64], gene_to_fam, means, &
                                          dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: above pi", 3)
        ! a gene without a family still needs an angle in range: there is no angle sentinel
        call compute_family_direction_rap(3, 1, [0.0_real64, 0.0_real64, -1.0_real64 - PI], gene_to_fam, means, &
                                          dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: out of range even without a family", 3)
        call compute_family_direction_rap(3, 1, [0.0_real64, NAN_VALUE, 0.0_real64], gene_to_fam, means, dispersions, counts, &
                                          status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "rap family: NaN angle", 3)
        call compute_family_direction_rap(3, 1, [0.0_real64, INF_VALUE, 0.0_real64], gene_to_fam, means, dispersions, counts, &
                                          status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "rap family: infinite angle", 3)
        call compute_family_direction_rap(3, 1, signed_angles, [1, 2, 0], means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: gene_to_fam above n_families", 4)
        call compute_family_direction_rap(3, 1, signed_angles, [1, -1, 0], means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: gene_to_fam below 1", 4)
        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          min_angular_dispersion=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: negative min", 9)
        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          max_angular_dispersion=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: negative max", 10)
        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          max_angular_dispersion=NAN_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "rap family: NaN max", 10)
        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          max_angular_dispersion=MAX_DISPERSION_LIMIT, ierr=ierr)
        call assert_err(ierr, ERR_OK, "rap family: max at the limit 5")
        call compute_family_direction_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap family: max above the limit 5", 10)
    end subroutine test_validation_compute_family_direction_rap

    subroutine test_validation_compute_angular_deviations_rap()
        real(real64) :: signed_angles(2), means(2), deviations(2)
        integer(int32) :: gene_to_fam(2), ierr

        signed_angles = [0.5_real64, PI]
        means = [PI, FAMILY_MEAN_ANGLE_SENTINEL]
        gene_to_fam = [1, 2]

        call compute_angular_deviations_rap(2, 2, signed_angles, means, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_OK, "rap deviations: valid input, the mean sentinel included")
        call compute_angular_deviations_rap(0, 2, signed_angles, means, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "rap deviations: n_genes 0", 1)
        call compute_angular_deviations_rap(2, 0, signed_angles, means, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "rap deviations: n_families 0", 2)
        call compute_angular_deviations_rap(2, 2, [-PI, 0.0_real64], means, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap deviations: angle -pi", 3)
        call compute_angular_deviations_rap(2, 2, [NAN_VALUE, 0.0_real64], means, gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "rap deviations: NaN angle", 3)
        call compute_angular_deviations_rap(2, 2, signed_angles, [-PI, 0.0_real64], gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap deviations: mean -pi", 4)
        call compute_angular_deviations_rap(2, 2, signed_angles, [3.5_real64, 0.0_real64], gene_to_fam, deviations, &
                                            ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap deviations: mean above pi", 4)
        call compute_angular_deviations_rap(2, 2, signed_angles, [-1.0_real64, INF_VALUE], gene_to_fam, deviations, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "rap deviations: infinite mean", 4)
        call compute_angular_deviations_rap(2, 2, signed_angles, means, [3, 1], deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap deviations: gene_to_fam above n_families", 5)
        call compute_angular_deviations_rap(2, 2, signed_angles, means, [1, -1], deviations, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "rap deviations: gene_to_fam below 1", 5)
    end subroutine test_validation_compute_angular_deviations_rap

    subroutine test_validation_compute_relative_angular_deviations()
        real(real64) :: deviations(3), dispersions(2), relative(3)
        integer(int32) :: gene_to_fam(3), ierr

        deviations = [0.0_real64, PI, ANGULAR_DEVIATIONS_SENTINEL]
        dispersions = [0.0_real64, ANGULAR_DISPERSION_SENTINEL]
        gene_to_fam = [1, 2, 0]

        call compute_relative_angular_deviations(3, 2, deviations, dispersions, gene_to_fam, relative, ierr=ierr)
        call assert_err(ierr, ERR_OK, "relative: valid input, bounds and sentinels included")
        call compute_relative_angular_deviations(0, 2, deviations, dispersions, gene_to_fam, relative, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "relative: n_genes 0", 1)
        call compute_relative_angular_deviations(3, 0, deviations, dispersions, gene_to_fam, relative, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "relative: n_families 0", 2)
        call compute_relative_angular_deviations(3, 2, [0.0_real64, -0.5_real64, 0.0_real64], dispersions, gene_to_fam, &
                                                 relative, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "relative: negative deviation other than the sentinel", 3)
        call compute_relative_angular_deviations(3, 2, [0.0_real64, 3.2_real64, 0.0_real64], dispersions, gene_to_fam, &
                                                 relative, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "relative: deviation above pi", 3)
        call compute_relative_angular_deviations(3, 2, [0.0_real64, NAN_VALUE, 0.0_real64], dispersions, gene_to_fam, &
                                                 relative, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "relative: NaN deviation", 3)
        call compute_relative_angular_deviations(3, 2, deviations, [-0.5_real64, 1.0_real64], gene_to_fam, relative, &
                                                 ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "relative: negative dispersion other than the sentinel", 4)
        call compute_relative_angular_deviations(3, 2, deviations, [1.0_real64, INF_VALUE], gene_to_fam, relative, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "relative: infinite dispersion", 4)
        call compute_relative_angular_deviations(3, 2, deviations, dispersions, [1, 3, 0], relative, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "relative: gene_to_fam above n_families", 5)
        call compute_relative_angular_deviations(3, 2, deviations, dispersions, [1, -1, 0], relative, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "relative: gene_to_fam below 1", 5)
    end subroutine test_validation_compute_relative_angular_deviations

    subroutine test_validation_compute_angle_outlier_threshold()
        real(real64) :: values(3), threshold
        integer(int32) :: perm(3), ierr

        values = [0.0_real64, 2.0_real64, RELATIVE_SENTINEL]

        call compute_angle_outlier_threshold(3, values, threshold, quantile_level=1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "threshold: valid input, level 1 included")
        call compute_angle_outlier_threshold(0, values, threshold, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "threshold: n_genes 0", 1)
        call compute_angle_outlier_threshold(3, [0.0_real64, -0.5_real64, 1.0_real64], threshold, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "threshold: negative value other than the sentinel", 2)
        call compute_angle_outlier_threshold(3, [0.0_real64, NAN_VALUE, 1.0_real64], threshold, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "threshold: NaN value", 2)
        call compute_angle_outlier_threshold(3, [0.0_real64, INF_VALUE, 1.0_real64], threshold, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "threshold: infinite value", 2)
        call compute_angle_outlier_threshold(3, values, threshold, quantile_level=-0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "threshold: level below 0", 4)
        ! a percentage is not a level: 95 is rejected, not silently clamped
        call compute_angle_outlier_threshold(3, values, threshold, quantile_level=95.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "threshold: level 95 above 1", 4)
        call compute_angle_outlier_threshold(3, values, threshold, quantile_level=NAN_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "threshold: NaN level", 4)

        ! the expert tier numbers the work array too, so the level is argument 5
        call compute_angle_outlier_threshold_expert(3, values, perm, threshold, quantile_level=0.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "threshold expert: valid input, level 0 included")
        call compute_angle_outlier_threshold_expert(0, values, perm, threshold, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "threshold expert: n_genes 0", 1)
        call compute_angle_outlier_threshold_expert(3, [0.0_real64, -2.0_real64, 1.0_real64], perm, threshold, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "threshold expert: negative value", 2)
        call compute_angle_outlier_threshold_expert(3, values, perm, threshold, quantile_level=1.5_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "threshold expert: level above 1", 5)
    end subroutine test_validation_compute_angle_outlier_threshold

    subroutine test_validation_flag_angle_outliers()
        real(real64) :: values(2)
        logical(c_bool) :: is_outlier(2)
        integer(int32) :: ierr

        values = [0.0_real64, RELATIVE_SENTINEL]

        call flag_angle_outliers(2, values, -1.0_real64, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_OK, "flag: valid input, any finite threshold")
        call flag_angle_outliers(0, values, 1.0_real64, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "flag: n_genes 0", 1)
        call flag_angle_outliers(2, [-0.5_real64, 0.0_real64], 1.0_real64, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "flag: negative value other than the sentinel", 2)
        call flag_angle_outliers(2, [NAN_VALUE, 0.0_real64], 1.0_real64, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "flag: NaN value", 2)
        call flag_angle_outliers(2, values, NAN_VALUE, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "flag: NaN threshold", 3)
        call flag_angle_outliers(2, values, INF_VALUE, is_outlier, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "flag: infinite threshold", 3)
    end subroutine test_validation_flag_angle_outliers

    subroutine test_validation_detect_angle_outliers()
        real(real64) :: expression_vectors(2, 3), bad(2, 3), directions(2, 1), dispersions(1), relative(3), threshold
        real(real64) :: tmp_deviations(3)
        integer(int32) :: gene_to_fam(3), counts(1), status(1), gene_status(3), tmp_perm(3), ierr
        logical(c_bool) :: is_outlier(3)

        expression_vectors = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64], [2, 3])
        gene_to_fam = [1, 1, 1]

        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect: valid input")
        call detect_angle_outliers(0, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "detect: n_axes 0", 1)
        call detect_angle_outliers(2, 0, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "detect: n_genes 0", 2)
        call detect_angle_outliers(2, 3, 0, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "detect: n_families 0", 3)
        bad = expression_vectors
        bad(1, 3) = NAN_VALUE
        call detect_angle_outliers(2, 3, 1, bad, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect: NaN expression", 4)
        bad(1, 3) = INF_VALUE
        call detect_angle_outliers(2, 3, 1, bad, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect: infinite expression", 4)
        call detect_angle_outliers(2, 3, 1, expression_vectors, [1, 2, 1], directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: gene_to_fam above n_families", 5)
        call detect_angle_outliers(2, 3, 1, expression_vectors, [1, -1, 1], directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: gene_to_fam below 1", 5)
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, quantile_level=95.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: level 95 above 1", 14)
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, quantile_level=-0.5_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: level below 0", 14)
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, min_angular_dispersion=-1.0_real64, &
                                   ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: negative min", 15)
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, max_angular_dispersion=NAN_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect: NaN max", 16)
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, max_angular_dispersion=MAX_DISPERSION_LIMIT, &
                                   ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect: max at the limit 5")
        call detect_angle_outliers(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, &
                                   max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect: max above the limit 5", 16)

        ! the expert tier numbers its two work arrays, so the optionals move by two
        call detect_angle_outliers_expert(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect expert: valid input")
        call detect_angle_outliers_expert(2, 3, 1, bad, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect expert: infinite expression", 4)
        call detect_angle_outliers_expert(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          quantile_level=2.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect expert: level above 1", 16)
        call detect_angle_outliers_expert(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          min_angular_dispersion=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect expert: negative min", 17)
        call detect_angle_outliers_expert(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          max_angular_dispersion=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect expert: negative max", 18)
        call detect_angle_outliers_expert(2, 3, 1, expression_vectors, gene_to_fam, directions, dispersions, counts, &
                                          status, relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                          max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect expert: max above the limit 5", 18)
    end subroutine test_validation_detect_angle_outliers

    subroutine test_validation_detect_angle_outliers_rap()
        real(real64) :: signed_angles(3), means(1), dispersions(1), relative(3), threshold, tmp_deviations(3)
        integer(int32) :: gene_to_fam(3), counts(1), status(1), gene_status(3), tmp_perm(3), ierr
        logical(c_bool) :: is_outlier(3)

        signed_angles = [0.1_real64, 0.2_real64, PI]
        gene_to_fam = [1, 1, 1]

        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect rap: valid input")
        call detect_angle_outliers_rap(0, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "detect rap: n_genes 0", 1)
        call detect_angle_outliers_rap(3, 0, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "detect rap: n_families 0", 2)
        call detect_angle_outliers_rap(3, 1, [0.1_real64, -PI, 0.3_real64], gene_to_fam, means, dispersions, counts, &
                                       status, relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: angle -pi", 3)
        call detect_angle_outliers_rap(3, 1, [0.1_real64, 4.0_real64, 0.3_real64], gene_to_fam, means, dispersions, &
                                       counts, status, relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: angle above pi", 3)
        call detect_angle_outliers_rap(3, 1, [0.1_real64, NAN_VALUE, 0.3_real64], gene_to_fam, means, dispersions, counts, &
                                       status, relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect rap: NaN angle", 3)
        call detect_angle_outliers_rap(3, 1, signed_angles, [1, 2, 1], means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: gene_to_fam above n_families", 4)
        call detect_angle_outliers_rap(3, 1, signed_angles, [1, -5, 1], means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: gene_to_fam below 1", 4)
        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, quantile_level=95.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: level 95 above 1", 13)
        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, min_angular_dispersion=NAN_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect rap: NaN min", 14)
        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, max_angular_dispersion=-0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: negative max", 15)
        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, max_angular_dispersion=MAX_DISPERSION_LIMIT, &
                                       ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect rap: max at the limit 5")
        call detect_angle_outliers_rap(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, &
                                       max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap: max above the limit 5", 15)

        call detect_angle_outliers_rap_expert(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                              relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                              ierr=ierr)
        call assert_err(ierr, ERR_OK, "detect rap expert: valid input")
        call detect_angle_outliers_rap_expert(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                              relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                              quantile_level=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap expert: level below 0", 15)
        call detect_angle_outliers_rap_expert(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                              relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                              min_angular_dispersion=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap expert: negative min", 16)
        call detect_angle_outliers_rap_expert(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                              relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                              max_angular_dispersion=INF_VALUE, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "detect rap expert: infinite max", 17)
        call detect_angle_outliers_rap_expert(3, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                              relative, threshold, is_outlier, gene_status, tmp_deviations, tmp_perm, &
                                              max_angular_dispersion=nearest(MAX_DISPERSION_LIMIT, 1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "detect rap expert: max above the limit 5", 17)
    end subroutine test_validation_detect_angle_outliers_rap

    !> The status codes are distinct from each other, none is zero -- the "nothing to report"
    !| status -- and none carries an argument position.
    subroutine test_status_codes_are_distinct_statuses()
        integer(int32), parameter :: codes(5) = [STAT_NO_STABLE_DIRECTION, STAT_NO_ANGULAR_VARIATION, &
                                                 STAT_TOO_FEW_MEMBERS, STAT_NO_FAMILY, STAT_ZERO_VECTOR]
        integer(int32) :: i_code, j_code

        do i_code = 1, size(codes)
            call assert_not_equal_int(codes(i_code), ERR_OK, "status code is not zero")
            call assert_equal_int(get_err_code(codes(i_code)), codes(i_code), "status code carries no argument position")
            do j_code = i_code + 1, size(codes)
                call assert_not_equal_int(codes(i_code), codes(j_code), "status codes distinct")
            end do
        end do
    end subroutine test_status_codes_are_distinct_statuses

    ! ------------------------------------------------------------------------------------------
    ! the noise floor, and the decisions it replaces
    ! ------------------------------------------------------------------------------------------

    !> Sines and cosines that cancel to exactly zero: no mean angle, not atan2(0, 0) = 0.
    !|
    !| The angles x, -x, y, -y with y = pi - x have sines summing to exactly 0 (the library sine is
    !| odd) and cosines summing to 2 cos x + 2 cos y, which is exactly 0 wherever the library gives
    !| cos y = -cos x to the bit. Such an x is searched for, since which ones qualify depends on
    !| the library; that at least one qualifies is asserted too.
    subroutine test_rap_exact_cancellation()
        real(real64) :: signed_angles(4), family_mean_angles(1), angular_dispersions(1), x, y
        integer(int32) :: gene_to_fam(4), member_counts(1), status(1), ierr, i_candidate
        logical :: found

        found = .false.
        do i_candidate = 1, 150
            x = 0.01_real64*real(i_candidate, real64)
            y = PI - x
            ! summed in the order the implementation sums them; `.not. abs(v) > 0` is `v == 0`
            if (.not. (abs(sin(-x) + sin(x)) > 0.0_real64 .or. abs(sin(-y) + sin(y)) > 0.0_real64 &
                       .or. abs((cos(x) + cos(-x)) + cos(y) + cos(-y)) > 0.0_real64)) then
                found = .true.
                exit
            end if
        end do
        call assert_true(found, "rap exact cancellation: an angle whose cosines cancel exactly exists")
        if (.not. found) return

        signed_angles(1) = x
        signed_angles(2) = -x
        signed_angles(3) = y
        signed_angles(4) = -y
        gene_to_fam = 1
        family_mean_angles = UNWRITTEN
        angular_dispersions = UNWRITTEN
        call compute_family_direction_rap(4, 1, signed_angles, gene_to_fam, family_mean_angles, angular_dispersions, &
                                          member_counts, status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "rap exact cancellation: ierr")
        call assert_equal_int(status(1), STAT_NO_STABLE_DIRECTION, "rap exact cancellation: no stable direction")
        call assert_equal_real(family_mean_angles(1), FAMILY_MEAN_ANGLE_SENTINEL, 0.0_real64, &
                               "rap exact cancellation: no mean angle")
        call assert_equal_real(angular_dispersions(1), ANGULAR_DISPERSION_SENTINEL, 0.0_real64, &
                               "rap exact cancellation: no dispersion")
    end subroutine test_rap_exact_cancellation

    !> Identical and proportional members stay within the noise floor at any size: 10000 axes,
    !| 10000 members, 10000 signed angles at a generic angle -- no angular variation, with the
    !| default minimum of 0.
    subroutine test_noise_floor_large_families()
        integer(int32), parameter :: n_wide = 10000, n_many = 10000
        real(real64), allocatable :: wide(:, :), many(:, :), wide_directions(:, :), many_directions(:, :)
        real(real64), allocatable :: angles(:)
        integer(int32), allocatable :: gene_to_fam_many(:), gene_to_fam_angles(:)
        real(real64) :: dispersions(1), means(1), scale
        integer(int32) :: counts(1), status(1), ierr, i_axis, i_gene

        ! 3 members of 10000 axes: the same generic vector at scales 1, 1e200 and 3e-150
        allocate (wide(n_wide, 3), wide_directions(n_wide, 1))
        do i_axis = 1, n_wide
            wide(i_axis, 1) = 1.0_real64 + real(mod(37*i_axis, 101), real64)/7.0_real64
        end do
        wide(:, 2) = 1.0e200_real64*wide(:, 1)
        wide(:, 3) = 3.0e-150_real64*wide(:, 1)
        call compute_family_direction(n_wide, 3, 1, wide, [1, 1, 1], wide_directions, dispersions, counts, status, &
                                      ierr=ierr)
        call assert_err(ierr, ERR_OK, "noise floor, 10000 axes: ierr")
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "noise floor, 10000 axes: no angular variation")
        call assert_all_finite(wide_directions, n_wide, "noise floor, 10000 axes: direction")

        ! 10000 members of 3 axes: identical, then proportional over 1e-300..1e300
        allocate (many(3, n_many), many_directions(3, 1), gene_to_fam_many(n_many))
        gene_to_fam_many = 1
        many(1, :) = 0.3_real64
        many(2, :) = -1.7_real64
        many(3, :) = 2.9_real64
        call compute_family_direction(3, n_many, 1, many, gene_to_fam_many, many_directions, dispersions, counts, &
                                      status, ierr=ierr)
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "noise floor, 10000 identical: no variation")
        do i_gene = 1, n_many
            scale = 10.0_real64**(real(mod(97*i_gene, 601) - 300, real64))
            many(1, i_gene) = scale*0.3_real64
            many(2, i_gene) = scale*(-1.7_real64)
            many(3, i_gene) = scale*2.9_real64
        end do
        call compute_family_direction(3, n_many, 1, many, gene_to_fam_many, many_directions, dispersions, counts, &
                                      status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "noise floor, 10000 proportional: ierr")
        call assert_equal_int(counts(1), n_many, "noise floor, 10000 proportional: all members")
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "noise floor, 10000 proportional: no variation")

        ! 10000 signed angles of 0.7
        allocate (angles(n_many), gene_to_fam_angles(n_many))
        angles = 0.7_real64
        gene_to_fam_angles = 1
        call compute_family_direction_rap(n_many, 1, angles, gene_to_fam_angles, means, dispersions, counts, status, &
                                          ierr=ierr)
        call assert_err(ierr, ERR_OK, "noise floor, 10000 angles: ierr")
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "noise floor, 10000 angles: no variation")
        ! the sums of 10000 sines and cosines turn the mean by up to (m - 1) eps
        call assert_equal_real(means(1), 0.7_real64, real(n_many, real64)*epsilon(1.0_real64), &
                               "noise floor, 10000 angles: mean")
    end subroutine test_noise_floor_large_families

    !> pi and -pi + ulp are 4.4e-16 rad apart around the circle: no angular variation, and the
    !| mean is the direction pi -- atan2 of these sums returns exactly -pi, which must be wrapped.
    subroutine test_rap_pi_mixed_with_minus_pi()
        real(real64) :: signed_angles(4), means(1), dispersions(1)
        integer(int32) :: gene_to_fam(4), counts(1), status(1), ierr

        signed_angles = [PI, nearest(-PI, 1.0_real64), PI, nearest(-PI, 1.0_real64)]
        gene_to_fam = 1
        call compute_family_direction_rap(4, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "pi and -pi + ulp: ierr")
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "pi and -pi + ulp: no angular variation")
        call assert_true(means(1) > -PI .and. means(1) <= PI, "pi and -pi + ulp: mean inside (-pi, pi]")
        call assert_equal_real(abs(wrap_angle(means(1) - PI)), 0.0_real64, TOL_ULP, "pi and -pi + ulp: mean is pi")
    end subroutine test_rap_pi_mixed_with_minus_pi

    !> Deviations {-d, 0, d}: 1 - R = (2/3)(1 - cos d) = (4/3) sin^2(d/2), so
    !| sigma = sqrt(-2 ln(1 - (4/3) sin^2(d/2))) = sqrt(8/3) sin(d/2) (1 + O(d^2)).
    !|
    !| d = 1e-10 and 1e-13 must come out to 1e-9 relative. Both lie above the noise floor of three
    !| members in two axes, (3 + 3 + 16) eps = 4.9e-15, by 17x at 1e-13; each would vanish if
    !| 1 - cos d were taken instead of 2 sin^2(d/2), or ln(1 - x) without its series (both round
    !| to 0 here), and 1e-13 falls below a floor a hundred times larger. The inputs make the
    !| deviations exact: the members are (1, -d), (1, 0), (1, d) -- cos d rounds to 1 -- whose sum
    !| (3, 0) is exact, and the angles -d, 0, d, whose sines cancel exactly. d = 1e-17 is below any
    !| rounding floor and counts as no angular variation.
    subroutine test_tiny_dispersion_above_and_below_floor()
        real(real64), parameter :: accepted(2) = [1.0e-10_real64, 1.0e-13_real64], d_rejected = 1.0e-17_real64
        real(real64) :: vectors(2, 3), directions(2, 1), dispersions(1), means(1), angles(3), d, sigma
        integer(int32) :: counts(1), status(1), ierr, i_d

        do i_d = 1, size(accepted)
            d = accepted(i_d)
            sigma = sqrt(8.0_real64/3.0_real64)*sin(0.5_real64*d)

            vectors(1, :) = 1.0_real64
            vectors(2, 1) = -d
            vectors(2, 2) = 0.0_real64
            vectors(2, 3) = d
            call compute_family_direction(2, 3, 1, vectors, [1, 1, 1], directions, dispersions, counts, status, &
                                          ierr=ierr)
            call assert_err(ierr, ERR_OK, "tiny dispersion: ierr")
            call assert_equal_int(status(1), ERR_OK, "tiny dispersion: accepted at d = "//d)
            call assert_equal_real(dispersions(1)/sigma, 1.0_real64, 1.0e-9_real64, "tiny dispersion: sigma at d = "//d)

            angles(1) = -d
            angles(2) = 0.0_real64
            angles(3) = d
            call compute_family_direction_rap(3, 1, angles, [1, 1, 1], means, dispersions, counts, status, ierr=ierr)
            call assert_equal_int(status(1), ERR_OK, "rap tiny dispersion: accepted at d = "//d)
            call assert_equal_real(dispersions(1)/sigma, 1.0_real64, 1.0e-9_real64, "rap tiny dispersion: sigma at d = "//d)
        end do

        vectors(1, :) = 1.0_real64
        vectors(2, :) = [-d_rejected, 0.0_real64, d_rejected]
        call compute_family_direction(2, 3, 1, vectors, [1, 1, 1], directions, dispersions, counts, status, ierr=ierr)
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "tiny dispersion 1e-17: below the floor")
        angles = [-d_rejected, 0.0_real64, d_rejected]
        call compute_family_direction_rap(3, 1, angles, [1, 1, 1], means, dispersions, counts, status, ierr=ierr)
        call assert_equal_int(status(1), STAT_NO_ANGULAR_VARIATION, "rap tiny dispersion 1e-17: below the floor")
    end subroutine test_tiny_dispersion_above_and_below_floor

    !> A nearly cancelling resultant is no stable direction, and like an exactly cancelling one it
    !| has no direction: RAP {0, pi/2, pi, -pi/2} (the sums are rounding residue, their atan2 an
    !| arbitrary angle) and spherical {e1, -e1, 2 e2, -3 e2} (unit vectors e1, -e1, e2, -e2, which
    !| a fast floating-point model need not cancel exactly). Their genes then get no deviation.
    subroutine test_nearly_cancelling_resultant()
        real(real64) :: signed_angles(4), vectors(2, 4), means(1), directions(2, 1), dispersions(1)
        real(real64) :: deviations(4), relative(4), threshold
        integer(int32) :: gene_to_fam(4), counts(1), status(1), gene_status(4), ierr
        logical(c_bool) :: is_outlier(4)

        gene_to_fam = 1
        signed_angles = [0.0_real64, PI/2, PI, -PI/2]
        means = UNWRITTEN
        call compute_family_direction_rap(4, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          ierr=ierr)
        call assert_err(ierr, ERR_OK, "nearly cancelling rap: ierr")
        call assert_equal_int(status(1), STAT_NO_STABLE_DIRECTION, "nearly cancelling rap: no stable direction")
        call assert_equal_real(means(1), FAMILY_MEAN_ANGLE_SENTINEL, 0.0_real64, "nearly cancelling rap: no mean angle")
        call compute_angular_deviations_rap(4, 1, signed_angles, means, gene_to_fam, deviations, ierr=ierr)
        call assert_equal_array_real(deviations, spread(ANGULAR_DEVIATIONS_SENTINEL, 1, 4), 4, 0.0_real64, &
                                     "nearly cancelling rap: no deviations")
        call detect_angle_outliers_rap(4, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, relative, &
                                       threshold, is_outlier, gene_status, ierr=ierr)
        call assert_true(all(gene_status == STAT_NO_STABLE_DIRECTION), "nearly cancelling rap: gene status")
        call assert_false(any(is_outlier), "nearly cancelling rap: nothing flagged")

        vectors = reshape([1.0_real64, 0.0_real64, -1.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, &
                           0.0_real64, -3.0_real64], [2, 4])
        directions = UNWRITTEN
        call compute_family_direction(2, 4, 1, vectors, gene_to_fam, directions, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "nearly cancelling spherical: ierr")
        call assert_equal_int(status(1), STAT_NO_STABLE_DIRECTION, "nearly cancelling spherical: no stable direction")
        call assert_equal_array_real(directions(:, 1), [0.0_real64, 0.0_real64], 2, 0.0_real64, &
                                     "nearly cancelling spherical: zero vector")
        call compute_angular_deviations(2, 4, 1, vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_equal_array_real(deviations, spread(ANGULAR_DEVIATIONS_SENTINEL, 1, 4), 4, 0.0_real64, &
                                     "nearly cancelling spherical: no deviations")
    end subroutine test_nearly_cancelling_resultant

    !> A vector whose components are all subnormal counts as a zero vector, whether the build
    !| flushes subnormals to zero or not: three copies of (0, 0, -1e-310, 0, 0, 0) used to give NaN
    !| with ifx's fast floating-point model and gradual underflow (as under Python or R). Where
    !| the underflow mode can be switched, both modes are tried, and the mode is restored after.
    subroutine test_subnormal_vector_is_zero_vector()
        logical :: gradual_at_start, gradual
        integer :: i_mode

        call ieee_get_underflow_mode(gradual_at_start)
        do i_mode = 1, 2
            gradual = i_mode == 1
            if (ieee_support_underflow_control(1.0_real64)) then
                call ieee_set_underflow_mode(gradual)
            else if (i_mode == 2) then
                exit  ! the mode cannot be switched: the build's own mode was tried once
            end if
            call check_subnormal_vector(gradual)
        end do
        if (ieee_support_underflow_control(1.0_real64)) call ieee_set_underflow_mode(gradual_at_start)
    end subroutine test_subnormal_vector_is_zero_vector

    subroutine check_subnormal_vector(gradual)
        logical, intent(in) :: gradual
        character(len=*), parameter :: SUBNORMAL_LABEL = "subnormal vector: "
        real(real64) :: vectors(6, 4), directions(6, 1), dispersions(1), deviations(4), relative(4), threshold
        integer(int32) :: gene_to_fam(4), counts(1), status(1), gene_status(4), ierr
        logical(c_bool) :: is_outlier(4)
        character(len=16) :: mode

        mode = merge("gradual         ", "flush to zero   ", gradual)
        vectors = 0.0_real64
        vectors(3, 1:3) = -1.0e-310_real64
        vectors(:, 4) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64]
        gene_to_fam = 1

        call compute_family_direction(6, 4, 1, vectors, gene_to_fam, directions, dispersions, counts, status, ierr=ierr)
        call assert_err(ierr, ERR_OK, SUBNORMAL_LABEL//trim(mode)//": ierr")
        call assert_equal_int(counts(1), 1, SUBNORMAL_LABEL//trim(mode)//": not counted as members")
        call assert_equal_int(status(1), STAT_TOO_FEW_MEMBERS, SUBNORMAL_LABEL//trim(mode)//": too few members")
        call assert_all_finite(directions, 6, SUBNORMAL_LABEL//trim(mode)//": direction")
        call assert_all_finite(dispersions, 1, SUBNORMAL_LABEL//trim(mode)//": dispersion")

        ! against a real direction, a subnormal vector has no angle
        directions(:, 1) = [1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        call compute_angular_deviations(6, 4, 1, vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call assert_equal_array_real(deviations(1:3), spread(ANGULAR_DEVIATIONS_SENTINEL, 1, 3), 3, 0.0_real64, &
                                     SUBNORMAL_LABEL//trim(mode)//": no angle")
        call assert_all_finite(deviations, 4, SUBNORMAL_LABEL//trim(mode)//": angles")

        call detect_angle_outliers(6, 4, 1, vectors, gene_to_fam, directions, dispersions, counts, status, relative, &
                                   threshold, is_outlier, gene_status, ierr=ierr)
        call assert_equal_array_int(gene_status(1:3), [STAT_ZERO_VECTOR, STAT_ZERO_VECTOR, STAT_ZERO_VECTOR], 3, &
                                    SUBNORMAL_LABEL//trim(mode)//": gene status")
        call assert_all_finite(relative, 4, SUBNORMAL_LABEL//trim(mode)//": relative deviations")
    end subroutine check_subnormal_vector

    !> m - 1 identical members and one at angle a: the direction turns by phi = atan2(sin a,
    !| m - 1 + cos a) towards it, R = ((m - 1) cos phi + cos(a - phi))/m, and the deviating
    !| member's relative angular deviation (a - phi)/sigma is about sqrt(m).
    subroutine test_single_deviating_member()
        integer(int32), parameter :: m = 100
        real(real64), parameter :: a = 0.3_real64
        real(real64) :: vectors(2, m), directions(2, 1), dispersions(1), deviations(m), relative(m)
        real(real64) :: phi, sigma, expected_relative, threshold
        integer(int32) :: gene_to_fam(m), counts(1), status(1), ierr
        logical(c_bool) :: is_outlier(m)

        vectors(1, :) = 1.0_real64
        vectors(2, :) = 0.0_real64
        vectors(:, m) = [cos(a), sin(a)]
        gene_to_fam = 1
        phi = atan2(sin(a), real(m - 1, real64) + cos(a))
        sigma = sqrt(-2.0_real64*log((real(m - 1, real64)*cos(phi) + cos(a - phi))/real(m, real64)))
        expected_relative = (a - phi)/sigma

        call compute_family_direction(2, m, 1, vectors, gene_to_fam, directions, dispersions, counts, status, ierr=ierr)
        call assert_equal_int(status(1), ERR_OK, "single deviating member: family accepted")
        call assert_equal_real(dispersions(1)/sigma, 1.0_real64, 1.0e-9_real64, "single deviating member: sigma")
        call compute_angular_deviations(2, m, 1, vectors, directions, gene_to_fam, deviations, ierr=ierr)
        call compute_relative_angular_deviations(m, 1, deviations, dispersions, gene_to_fam, relative, ierr=ierr)
        call assert_equal_real(relative(m)/expected_relative, 1.0_real64, 1.0e-9_real64, &
                               "single deviating member: relative deviation (a - phi)/sigma")
        call assert_in_range_real(relative(m)/sqrt(real(m, real64)), 0.99_real64, 1.01_real64, &
                                  "single deviating member: about sqrt(m)")
        call compute_angle_outlier_threshold(m, relative, threshold, ierr=ierr)
        call flag_angle_outliers(m, relative, threshold, is_outlier, ierr=ierr)
        ! the 99 identical members tie at the threshold and are not flagged; only the deviating one is
        call assert_true(is_outlier(m), "single deviating member: flagged")
        call assert_equal_int(count(is_outlier), 1, "single deviating member: the only one flagged")
    end subroutine test_single_deviating_member

    !> With the default minimum (0), a family of identical members still has no angular
    !| variation: its genes report it and none is flagged, while a usable family beside it is ranked.
    subroutine test_detect_identical_family_default_min()
        integer(int32), parameter :: n_genes = 6
        real(real64) :: vectors(2, n_genes), directions(2, 2), dispersions(2), relative(n_genes), threshold
        integer(int32) :: gene_to_fam(n_genes), counts(2), status(2), gene_status(n_genes), ierr
        logical(c_bool) :: is_outlier(n_genes)

        ! family 1: (2, 5) three times; family 2: e1, (1, 1), e2 (see test_detect_angle_outliers)
        vectors = reshape([2.0_real64, 5.0_real64, 2.0_real64, 5.0_real64, 2.0_real64, 5.0_real64, &
                           1.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64, 1.0_real64], [2, n_genes])
        gene_to_fam = [1, 1, 1, 2, 2, 2]
        is_outlier = .true.
        call detect_angle_outliers(2, n_genes, 2, vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "identical family, default min: ierr")
        call assert_equal_array_int(status, [STAT_NO_ANGULAR_VARIATION, ERR_OK], 2, "identical family, default min: status")
        call assert_equal_array_int(gene_status(1:3), [STAT_NO_ANGULAR_VARIATION, STAT_NO_ANGULAR_VARIATION, &
                                    STAT_NO_ANGULAR_VARIATION], 3, "identical family, default min: gene status")
        call assert_false(any(is_outlier(1:3)), "identical family, default min: none flagged")
        call assert_equal_array_int(gene_status(4:6), [ERR_OK, ERR_OK, ERR_OK], 3, &
                                    "identical family, default min: the other family ranked")
    end subroutine test_detect_identical_family_default_min

    !> 999 identical genes and one 0.3 rad off: sigma is 9.45e-3, which a minimum of 0.01 used
    !| to reject, so the clear outlier was never flagged. With the default minimum of 0 it is.
    subroutine test_detect_999_plus_1_not_masked()
        integer(int32), parameter :: n_genes = 1000
        real(real64), parameter :: a = 0.3_real64
        real(real64) :: vectors(2, n_genes), directions(2, 1), dispersions(1), relative(n_genes), threshold
        real(real64) :: phi, sigma
        integer(int32) :: gene_to_fam(n_genes), counts(1), status(1), gene_status(n_genes), ierr
        logical(c_bool) :: is_outlier(n_genes)

        vectors(1, :) = 1.0_real64
        vectors(2, :) = 0.0_real64
        vectors(:, n_genes) = [cos(a), sin(a)]
        gene_to_fam = 1
        ! see test_single_deviating_member: sigma = 9.45e-3
        phi = atan2(sin(a), real(n_genes - 1, real64) + cos(a))
        sigma = sqrt(-2.0_real64*log((real(n_genes - 1, real64)*cos(phi) + cos(a - phi))/real(n_genes, real64)))

        call detect_angle_outliers(2, n_genes, 1, vectors, gene_to_fam, directions, dispersions, counts, status, &
                                   relative, threshold, is_outlier, gene_status, ierr=ierr)
        call assert_err(ierr, ERR_OK, "999 + 1: ierr")
        call assert_equal_int(status(1), ERR_OK, "999 + 1: family accepted by the default minimum")
        call assert_in_range_real(sigma, 9.4e-3_real64, 9.5e-3_real64, "999 + 1: the hand-derived sigma is 9.45e-3")
        call assert_equal_real(dispersions(1)/sigma, 1.0_real64, 1.0e-9_real64, "999 + 1: sigma")
        ! the 999 identical genes tie at the threshold, so only the outlier lies above it
        call assert_true(is_outlier(n_genes), "999 + 1: the outlier is flagged")
        call assert_equal_int(count(is_outlier), 1, "999 + 1: only the outlier is flagged")
    end subroutine test_detect_999_plus_1_not_masked

    !> The documented defaults are what an absent argument means.
    subroutine test_published_defaults()
        real(real64) :: vectors(2, 7), directions(2, 2), dispersions(2), directions_explicit(2, 2)
        real(real64) :: dispersions_explicit(2), threshold, threshold_explicit, values(4)
        integer(int32) :: counts(2), status(2), status_explicit(2), ierr

        ! family 1 has sigma sqrt(ln 2), family 2 sigma sqrt(2 ln 3) (see test_family_direction_exact)
        vectors = reshape([1.0_real64, 0.0_real64, 3.0_real64, 0.0_real64, 0.0_real64, 2.0_real64, 0.0_real64, 0.5_real64, &
                           1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, -1.0_real64, 0.0_real64], [2, 7])
        call compute_family_direction(2, 7, 2, vectors, [1, 1, 1, 1, 2, 2, 2], directions, dispersions, counts, status, &
                                      ierr=ierr)
        call compute_family_direction(2, 7, 2, vectors, [1, 1, 1, 1, 2, 2, 2], directions_explicit, &
                                      dispersions_explicit, counts, status_explicit, &
                                      min_angular_dispersion=MIN_DISPERSION_DEFAULT, &
                                      max_angular_dispersion=MAX_DISPERSION_DEFAULT, ierr=ierr)
        call assert_equal_array_int(status, status_explicit, 2, "defaults: dispersion bounds")
        call assert_equal_array_real(dispersions, dispersions_explicit, 2, 0.0_real64, "defaults: dispersions")

        values = [0.5_real64, 4.0_real64, 1.0_real64, 2.0_real64]
        call compute_angle_outlier_threshold(4, values, threshold, ierr=ierr)
        call compute_angle_outlier_threshold(4, values, threshold_explicit, quantile_level=0.95_real64, ierr=ierr)
        call assert_equal_real(threshold, threshold_explicit, 0.0_real64, "defaults: quantile level 0.95")
    end subroutine test_published_defaults

    !> Unit vectors that cancel out up to rounding leave a resultant of about 1e-16, whose
    !| dispersion sqrt(-2 ln R) is about 8.6 -- a "stable direction" for any maximum above that.
    !| The maximum is capped at 5 (R >= exp(-12.5)), so they are never one: RAP {t, t + pi/2,
    !| t + pi, t - pi/2} over 400 angles t, and spherical {u, -2u, v, -3v} over 400 pairs of
    !| generic directions, give STAT_NO_STABLE_DIRECTION at the default maximum and at the limit,
    !| and a maximum above the limit is refused.
    subroutine test_cancellation_never_stable_up_to_the_limit()
        integer(int32), parameter :: n_cases = 400
        real(real64) :: signed_angles(4), means(1), dispersions(1), vectors(3, 4), directions(3, 1), t
        real(real64) :: u(3), v(3), maxima(2)
        integer(int32) :: gene_to_fam(4), counts(1), status(1), ierr, i_case, i_max, n_stable

        gene_to_fam = 1
        maxima = [MAX_DISPERSION_DEFAULT, MAX_DISPERSION_LIMIT]
        do i_max = 1, size(maxima)
            n_stable = 0
            do i_case = 1, n_cases
                t = PI*(2.0_real64*real(i_case, real64)/real(n_cases, real64) - 1.0_real64) + 1.0e-3_real64
                signed_angles(1) = wrap_angle(t)
                signed_angles(2) = wrap_angle(t + PI/2)
                signed_angles(3) = wrap_angle(t + PI)
                signed_angles(4) = wrap_angle(t - PI/2)
                call compute_family_direction_rap(4, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                                  max_angular_dispersion=maxima(i_max), ierr=ierr)
                if (status(1) /= STAT_NO_STABLE_DIRECTION .or. ierr /= ERR_OK) n_stable = n_stable + 1
            end do
            call assert_equal_int(n_stable, 0, "cancelling rap families: none stable, max "//maxima(i_max))

            n_stable = 0
            do i_case = 1, n_cases
                u(1) = cos(0.37_real64*i_case)
                u(2) = sin(0.37_real64*i_case)
                u(3) = 0.3_real64 + 0.001_real64*i_case
                v(1) = 0.5_real64 - 0.002_real64*i_case
                v(2) = cos(1.3_real64*i_case)
                v(3) = sin(0.7_real64*i_case)
                vectors(:, 1) = u
                vectors(:, 2) = -2.0_real64*u
                vectors(:, 3) = v
                vectors(:, 4) = -3.0_real64*v
                call compute_family_direction(3, 4, 1, vectors, gene_to_fam, directions, dispersions, counts, status, &
                                              max_angular_dispersion=maxima(i_max), ierr=ierr)
                if (status(1) /= STAT_NO_STABLE_DIRECTION .or. ierr /= ERR_OK) n_stable = n_stable + 1
            end do
            call assert_equal_int(n_stable, 0, "cancelling spherical families: none stable, max "//maxima(i_max))
        end do

        ! the cap is what keeps them out: a maximum above it is refused
        call compute_family_direction_rap(4, 1, signed_angles, gene_to_fam, means, dispersions, counts, status, &
                                          max_angular_dispersion=100.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "cancelling families: max 100 refused", 10)
    end subroutine test_cancellation_never_stable_up_to_the_limit
end module mod_test_tox_get_outliers_by_angle
