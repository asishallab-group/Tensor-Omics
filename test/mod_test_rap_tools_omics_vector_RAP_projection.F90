!> Unit test suite for omics_vector_RAP_projection routine.
module mod_test_rap_tools_omics_vector_RAP_projection
    use asserts
    use tox_relative_axis_plane_tools
    use tox_relative_axis_plane_tools_impl, only: omics_vector_RAP_projection_impl
    use tox_normalization, only: normalize_unit_length
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_next_after
    use, intrinsic :: ieee_exceptions, only: ieee_status_type, ieee_get_status, ieee_set_status, &
                                             ieee_get_flag, ieee_set_flag, ieee_invalid
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_rap_tools_omics_vector_RAP_projection() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(12))

        all_tests(1) = test_case("test_omics_vector_RAP_projection_all_selected", test_all_selected)
        all_tests(2) = test_case("test_omics_vector_RAP_projection_one_axis_selected", test_one_axis_selected)
        all_tests(3) = test_case("test_omics_vector_RAP_projection_one_vector_selected", test_one_vector_selected)
        all_tests(4) = test_case("test_omics_vector_RAP_projection_constant_vector", test_constant_vector)
        all_tests(5) = test_case("test_omics_vector_RAP_projection_orthogonal_vector", test_orthogonal_vector)
        all_tests(6) = test_case("test_omics_vector_RAP_projection_no_axes", test_no_axes)
        all_tests(7) = test_case("test_omics_vector_RAP_projection_no_vectors", test_no_vectors)
        all_tests(8) = test_case("test_omics_vector_RAP_projection_mixed_selection", test_mixed_selection)
        all_tests(9) = test_case("test_omics_vector_RAP_projection_non_square_vecs", test_non_square_vecs)
        all_tests(10) = test_case("test_omics_vector_RAP_projection_equal_coordinates_give_exact_zeros", &
                                  test_equal_coordinates_give_exact_zeros)
        all_tests(11) = test_case("test_omics_vector_RAP_projection_equal_selected_coordinates_give_exact_zeros", &
                                  test_equal_selected_coordinates_give_exact_zeros)
        all_tests(12) = test_case("test_omics_vector_RAP_projection_nan_stays_quiet", test_nan_stays_quiet)
    end function get_all_tests_rap_tools_omics_vector_RAP_projection

    !> Wrapper function for the actual call of `call_omics_vector_RAP_projection`
    subroutine call_omics_vector_RAP_projection(test_name, vecs, axes_mask, vecs_mask, result_projections)

        character(len=*), intent(in) :: test_name
        real(real64), dimension(:, :), intent(in) :: vecs
        logical(c_bool), dimension(:), intent(in) :: axes_mask, vecs_mask
        real(real64), dimension(:, :), intent(out), optional :: result_projections

        integer(int32) :: n_axes, n_vecs, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: i_vec, ierr

        n_axes = size(axes_mask)
        n_vecs = size(vecs_mask)
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: "//trim(test_name))

        do i_vec = 1, n_selected_vecs
            call assert_equal_real( &
                sum(projections(:, i_vec))/real(n_selected_axes, real64), &
                0.0_real64, &
                1d-12, &
                "test_omics_vector_RAP_projection: Case '"//trim(test_name)//"' failed" &
                )
        end do

        if (present(result_projections)) then
            result_projections = projections
        end if
    end subroutine call_omics_vector_RAP_projection

    !> Test all axes and vectors are selected
    subroutine test_all_selected()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs, i_vec

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: all selected")

        ! Check the projection
        do i_vec = 1, n_selected_vecs
            call assert_equal_real( &
                sum(projections(:, i_vec))/real(n_selected_axes, real64), &
                0.0_real64, &
                1d-12, &
                "test_omics_vector_RAP_projection_all_selected: projection failed" &
                )
        end do
    end subroutine test_all_selected

    !> Test one axis and all vectors are selected
    subroutine test_one_axis_selected()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.false., .true., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: one axis selected")
    end subroutine test_one_axis_selected

    !> Test all axes and one vector are selected
    subroutine test_one_vector_selected()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .true., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: one vector selected")
    end subroutine test_one_vector_selected

    !> Test constant vector
    subroutine test_constant_vector()
        real(real64), dimension(3, 3) :: vecs
        real(real64), allocatable :: projections(:, :)
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64) :: expected_projection(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = 0.0_real64
        vecs(:, 1) = [5.0, 5.0, 5.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: constant vector")
        expected_projection = 0.0_real64
        call assert_equal_array_real( &
            projections(:, 1), &
            expected_projection, &
            3, &
            1d-12, &
            "test_omics_vector_RAP_projection_constant_vector: Expected zero vector" &
            )
    end subroutine test_constant_vector

    !> Test orthogonal vector
    subroutine test_orthogonal_vector()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = 0.0_real64
        vecs(:, 1) = [1.0, 0.0, -1.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: orthogonal vector")
    end subroutine test_orthogonal_vector

    !> Test no axes
    subroutine test_no_axes()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.false., .false., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, arg_pos=7_int32), "ierr should be ERR_EMPTY_INPUT for invalid input: no axes")
    end subroutine test_no_axes

    !> Test no vectors
    subroutine test_no_vectors()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .false., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, arg_pos=5_int32), "ierr should be ERR_EMPTY_INPUT for invalid input: no axes")
    end subroutine test_no_vectors

    !> Test mixed selection
    subroutine test_mixed_selection()
        real(real64), dimension(3, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        axes_mask = [.true., .false., .true.]
        vecs_mask = [.true., .false., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: mixed selection")
    end subroutine test_mixed_selection

    !> Test non-square vecs
    subroutine test_non_square_vecs()
        real(real64), dimension(4, 2) :: vecs
        logical(c_bool) :: axes_mask(4), vecs_mask(2)
        real(real64), allocatable :: projections(:, :)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        ! Define two 4D vectors
        vecs(:, 1) = [1.0, 2.0, 3.0, 4.0]
        vecs(:, 2) = [4.0, 3.0, 2.0, 1.0]

        ! Select all axes and both vectors
        axes_mask = [.true., .true., .true., .true.]
        vecs_mask = [.true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_vector_RAP_projection(vecs, 4, 2, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: non-square vecs")
    end subroutine test_non_square_vecs

    !> A vector whose coordinates are all equal projects to exact zeros, whatever its length or
    !| value, and normalize_unit_length then refuses it as the zero vector. Subtracting the rounded
    !| mean used to leave a residue for many of them ([0.7]*6 gives about -1.1e-16 per coordinate),
    !| which normalization turned into a unit vector off the plane.
    subroutine test_equal_coordinates_give_exact_zeros()
        real(real64) :: smallest_subnormal

        smallest_subnormal = ieee_next_after(0.0_real64, 1.0_real64)

        call check_equal_coordinates_give_exact_zeros("[0.7]*6", 0.7_real64, 6)
        call check_equal_coordinates_give_exact_zeros("[0.1]*3", 0.1_real64, 3)
        call check_equal_coordinates_give_exact_zeros("[0.1]*7", 0.1_real64, 7)
        call check_equal_coordinates_give_exact_zeros("[3.3]*13", 3.3_real64, 13)
        call check_equal_coordinates_give_exact_zeros("[1e300]*3", 1.0e300_real64, 3)
        call check_equal_coordinates_give_exact_zeros("[huge]*4", huge(1.0_real64), 4)
        call check_equal_coordinates_give_exact_zeros("[-2.5]*4", -2.5_real64, 4)
        call check_equal_coordinates_give_exact_zeros("[tiny]*5", tiny(1.0_real64), 5)
        call check_equal_coordinates_give_exact_zeros("[tiny/4]*3, subnormal", tiny(1.0_real64)/4, 3)
        call check_equal_coordinates_give_exact_zeros("[smallest subnormal]*7", smallest_subnormal, 7)
        call check_equal_coordinates_give_exact_zeros("[-smallest subnormal]*2", -smallest_subnormal, 2)
        call check_equal_coordinates_give_exact_zeros("[0]*3", 0.0_real64, 3)
        call check_equal_coordinates_give_exact_zeros("[5]*1, a single axis", 5.0_real64, 1)
    end subroutine test_equal_coordinates_give_exact_zeros

    !> Projects `[value]*n_axes` and checks the projection is exact zeros that normalization refuses.
    subroutine check_equal_coordinates_give_exact_zeros(case_name, value, n_axes)
        character(len=*), intent(in) :: case_name
        real(real64), intent(in) :: value
        integer(int32), intent(in) :: n_axes

        real(real64) :: vecs(n_axes, 1), projections(n_axes, 1), zeros(n_axes)
        logical(c_bool) :: axes_mask(n_axes), vecs_mask(1)
        integer(int32) :: ierr

        vecs = value
        axes_mask = .true.
        vecs_mask = .true.
        zeros = 0.0_real64
        projections = 1.0_real64

        call omics_vector_RAP_projection(vecs, n_axes, 1, vecs_mask, 1, axes_mask, n_axes, projections, ierr)
        call assert_err(ierr, ERR_OK, "equal coordinates "//case_name//": projection")
        call assert_identical_array_real(projections(:, 1), zeros, n_axes, &
                                         "equal coordinates "//case_name//": projection is not exact zeros")

        call normalize_unit_length(projections(:, 1), n_axes, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "equal coordinates "//case_name//": normalization must refuse")
    end subroutine check_equal_coordinates_give_exact_zeros

    !> Only the selected axes count: a vector equal on them projects to exact zeros whatever its
    !| other coordinates hold, beside vectors that are not and project as before.
    subroutine test_equal_selected_coordinates_give_exact_zeros()
        real(real64) :: vecs(4, 3), projections(3, 3), zeros(3)
        logical(c_bool) :: axes_mask(4), vecs_mask(3)
        integer(int32) :: ierr

        vecs(:, 1) = [0.7_real64, 9.0_real64, 0.7_real64, 0.7_real64]
        vecs(:, 2) = [1.0_real64, 2.0_real64, 4.0_real64, 1.0_real64]
        vecs(:, 3) = [0.1_real64, -3.0_real64, 0.1_real64, 0.1_real64]
        axes_mask = [.true., .false., .true., .true.]
        vecs_mask = .true.
        zeros = 0.0_real64

        call omics_vector_RAP_projection(vecs, 4, 3, vecs_mask, 3, axes_mask, 3, projections, ierr)
        call assert_err(ierr, ERR_OK, "equal selected coordinates: projection")
        call assert_identical_array_real(projections(:, 1), zeros, 3, "equal selected coordinates: first vector")
        call assert_identical_array_real(projections(:, 2), [-1.0_real64, 2.0_real64, -1.0_real64], 3, &
                                         "equal selected coordinates: the unequal vector")
        call assert_identical_array_real(projections(:, 3), zeros, 3, "equal selected coordinates: third vector")
    end subroutine test_equal_selected_coordinates_give_exact_zeros

    !> The comparison that finds equal coordinates must stay quiet on a NaN: no IEEE invalid flag,
    !| so a host halting on invalid is not stopped (the generated entry point refuses NaN first;
    !| this calls the implementation directly to reach the comparison with one).
    subroutine test_nan_stays_quiet()
        real(real64) :: vecs(3, 5), projections(3, 5), nan
        logical(c_bool) :: axes_mask(3), vecs_mask(5)
        type(ieee_status_type) :: status_on_entry
        logical :: invalid_signalled

        call ieee_get_status(status_on_entry)

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        vecs(:, 1) = [nan, nan, nan]
        vecs(:, 2) = [nan, 1.0_real64, 2.0_real64]
        vecs(:, 3) = [1.0_real64, 1.0_real64, nan]
        vecs(:, 4) = [1.0_real64, nan, 1.0_real64]
        vecs(:, 5) = [0.7_real64, 0.7_real64, 0.7_real64]
        axes_mask = .true.
        vecs_mask = .true.

        call ieee_set_flag(ieee_invalid, .false.)
        call omics_vector_RAP_projection_impl(vecs, 3, 5, vecs_mask, 5, axes_mask, 3, projections)
        call ieee_get_flag(ieee_invalid, invalid_signalled)
        call ieee_set_status(status_on_entry)

        call assert_false(invalid_signalled, "a NaN coordinate raised the IEEE invalid flag")
        call assert_identical_array_real(projections(:, 5), [0.0_real64, 0.0_real64, 0.0_real64], 3, &
                                         "nan stays quiet: the equal vector beside the NaNs")
    end subroutine test_nan_stays_quiet

end module mod_test_rap_tools_omics_vector_RAP_projection
