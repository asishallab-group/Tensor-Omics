!> Unit test suite for RAP Projection routines.
module mod_test_rap_tools_omics_field_RAP_projection
    use asserts
    use tox_relative_axis_plane_tools
    use tox_normalization, only: normalize_unit_length
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_next_after
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_rap_tools_omics_field_RAP_projection() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(11))

        all_tests(1) = test_case("test_omics_field_RAP_projection_all_selected", test_all_selected)
        all_tests(2) = test_case("test_omics_field_RAP_projection_one_axis_selected", test_one_axis_selected)
        all_tests(3) = test_case("test_omics_field_RAP_projection_one_vector_selected", test_one_vector_selected)
        all_tests(4) = test_case("test_omics_field_RAP_projection_constant_vector", test_constant_vector)
        all_tests(5) = test_case("test_omics_field_RAP_projection_orthogonal_vector", test_orthogonal_vector)
        all_tests(6) = test_case("test_omics_field_RAP_projection_no_axes", test_no_axes)
        all_tests(7) = test_case("test_omics_field_RAP_projection_no_vectors", test_no_vectors)
        all_tests(8) = test_case("test_omics_field_RAP_projection_mixed_selection", test_mixed_selection)
        all_tests(9) = test_case("test_omics_field_RAP_projection_non_square_vecs", test_non_square_vecs)
        all_tests(10) = test_case("test_omics_field_RAP_projection_concrete_example", test_concrete_example)
        all_tests(11) = test_case("test_omics_field_RAP_projection_equal_shifts_give_exact_zeros", &
                                  test_equal_shifts_give_exact_zeros)
    end function get_all_tests_rap_tools_omics_field_RAP_projection

    !> Test all axes and vectors are selected
    subroutine test_all_selected()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs, i_vec
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: all selected")

        do i_vec = 1, n_selected_vecs
            call assert_equal_real( &
                sum(projections(:, i_vec))/real(n_selected_axes, real64), &
                0.0_real64, &
                1d-12, &
                "test_omics_field_RAP_projection_all_selected: projection failed" &
                )
        end do
    end subroutine test_all_selected

    !> Test one axis and all vectors are selected
    subroutine test_one_axis_selected()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.false., .true., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: one axis selected")
    end subroutine test_one_axis_selected

    !> Test all axes and one vector are selected
    subroutine test_one_vector_selected()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .true., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: one vector selected")
    end subroutine test_one_vector_selected

    !> Test constant vector
    subroutine test_constant_vector()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        real(real64), allocatable :: projections(:, :)
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        real(real64) :: expected_projection(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = 0.0_real64
        vecs(:, 2, 1) = [5.0, 5.0, 5.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: constant vector")

        expected_projection = 0.0_real64
        call assert_equal_array_real( &
            projections(:, 1), &
            expected_projection, &
            3, &
            1d-12, &
            "test_omics_field_RAP_projection_constant_vector: Expected zero vector" &
            )
    end subroutine test_constant_vector

    !> Test orthogonal vector
    subroutine test_orthogonal_vector()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = 0.0_real64
        vecs(:, 1, 1) = [1.0, 0.0, -1.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: orthogonal vector")
    end subroutine test_orthogonal_vector

    !> Test no axes selected (error)
    subroutine test_no_axes()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.false., .false., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_true(ierr /= 0, "ierr should be nonzero for no axes selected")
    end subroutine test_no_axes

    !> Test no vectors selected (error)
    subroutine test_no_vectors()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .false., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_true(ierr /= 0, "ierr should be nonzero for no vectors selected")
    end subroutine test_no_vectors

    !> Test mixed selection
    subroutine test_mixed_selection()
        implicit none

        real(real64), dimension(3, 2, 3) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 2, 3])
        axes_mask = [.true., .false., .true.]
        vecs_mask = [.true., .false., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: mixed selection")
    end subroutine test_mixed_selection

    !> Test non-square vecs
    subroutine test_non_square_vecs()
        real(real64), dimension(4, 2, 1) :: vecs
        logical(c_bool) :: axes_mask(4), vecs_mask(1)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs(:, 1, 1) = [1.0, 2.0, 3.0, 4.0]
        vecs(:, 2, 1) = [4.0, 3.0, 2.0, 1.0]
        axes_mask = [.true., .true., .true., .true.]
        vecs_mask = [.true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 4, 1, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: non-square vecs")
    end subroutine test_non_square_vecs

    !> Test concrete example
    subroutine test_concrete_example()
        real(real64), dimension(3, 2, 1) :: vecs
        logical(c_bool) :: axes_mask(3), vecs_mask(1)
        real(real64), allocatable :: projections(:, :)
        real(real64), dimension(3) :: expected_projection_vec
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs(:, 1, 1) = [1.0_real64, -3.0_real64, 1.1_real64]
        vecs(:, 2, 1) = [3.0_real64, 6.0_real64, 2.2_real64]
        expected_projection_vec = vecs(:, 1, 1) - vecs(:, 2, 1)
        expected_projection_vec = expected_projection_vec - sum(expected_projection_vec)/3

        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection(vecs, 3, 1, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "ierr should be 0 for valid input: concrete example")

        call assert_equal_array_real( &
            projections(:, 1), &
            expected_projection_vec, &
            3, &
            1d-12, &
            "test_omics_field_RAP_projection_constant_vector: Calculated projection doesn't match expected" &
            )
    end subroutine test_concrete_example

    !> A field whose computed differences, origin minus target, are all equal projects to exact
    !| zeros, and normalize_unit_length then refuses it as the zero vector. Subtracting the rounded
    !| mean used to leave a residue for many of them, which normalization turned into a direction.
    subroutine test_equal_shifts_give_exact_zeros()
        real(real64) :: smallest_subnormal

        smallest_subnormal = ieee_next_after(0.0_real64, 1.0_real64)

        call check_equal_shifts_give_exact_zeros("0.7 - 0 on 6 axes", spread(0.7_real64, 1, 6), spread(0.0_real64, 1, 6))
        call check_equal_shifts_give_exact_zeros("0.1 - 0 on 7 axes", spread(0.1_real64, 1, 7), spread(0.0_real64, 1, 7))
        call check_equal_shifts_give_exact_zeros("3.3 - 0 on 13 axes", spread(3.3_real64, 1, 13), spread(0.0_real64, 1, 13))
        call check_equal_shifts_give_exact_zeros("1e300 - -1e300 on 3 axes", spread(1.0e300_real64, 1, 3), &
                                                 spread(-1.0e300_real64, 1, 3))
        call check_equal_shifts_give_exact_zeros("-2.5 - 0 on 4 axes", spread(-2.5_real64, 1, 4), spread(0.0_real64, 1, 4))
        call check_equal_shifts_give_exact_zeros("tiny - 0 on 5 axes", spread(tiny(1.0_real64), 1, 5), &
                                                 spread(0.0_real64, 1, 5))
        call check_equal_shifts_give_exact_zeros("smallest subnormal - 0 on 7 axes", spread(smallest_subnormal, 1, 7), &
                                                 spread(0.0_real64, 1, 7))
        ! distinct points whose differences all round to the same 0.7
        call check_equal_shifts_give_exact_zeros("distinct points, shift 0.7", [1.5_real64, 2.5_real64, 0.75_real64], &
                                                 [0.8_real64, 1.8_real64, 0.05_real64])
    end subroutine test_equal_shifts_give_exact_zeros

    !> Projects the field from `origin_point` to `target_point` and checks its projection is exact
    !| zeros that normalization refuses.
    subroutine check_equal_shifts_give_exact_zeros(case_name, origin_point, target_point)
        character(len=*), intent(in) :: case_name
        real(real64), intent(in) :: origin_point(:), target_point(:)

        real(real64) :: fields(size(origin_point), 2, 1), projections(size(origin_point), 1), zeros(size(origin_point))
        logical(c_bool) :: axes_mask(size(origin_point)), fields_mask(1)
        integer(int32) :: ierr, n_axes

        n_axes = size(origin_point)
        fields(:, 1, 1) = origin_point
        fields(:, 2, 1) = target_point
        axes_mask = .true.
        fields_mask = .true.
        zeros = 0.0_real64
        projections = 1.0_real64

        call omics_field_RAP_projection(fields, n_axes, 1, fields_mask, 1, axes_mask, n_axes, projections, ierr)
        call assert_err(ierr, ERR_OK, "equal shifts "//case_name//": projection")
        call assert_identical_array_real(projections(:, 1), zeros, n_axes, &
                                         "equal shifts "//case_name//": projection is not exact zeros")

        call normalize_unit_length(projections(:, 1), n_axes, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "equal shifts "//case_name//": normalization must refuse")
    end subroutine check_equal_shifts_give_exact_zeros

end module mod_test_rap_tools_omics_field_RAP_projection
