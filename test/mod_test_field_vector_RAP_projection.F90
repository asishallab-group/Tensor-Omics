! filepath: test/mod_test_field_vector_RAP_projection.f90
!> Unit test suite for RAP Projection routines.
module mod_test_rap_tools_omics_field_RAP_projection
    use asserts
    use tox_relative_axis_plane_tools
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_rap_tools_omics_field_RAP_projection() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))

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
    end function get_all_tests_rap_tools_omics_field_RAP_projection

    !> Test all axes and vectors are selected
    subroutine test_all_selected()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs, i_vec
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: all selected")

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

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.false., .true., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: one axis selected")
    end subroutine test_one_axis_selected

    !> Test all axes and one vector are selected
    subroutine test_one_vector_selected()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .true., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: one vector selected")
    end subroutine test_one_vector_selected

    !> Test constant vector
    subroutine test_constant_vector()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        real(real64), allocatable :: projections(:, :)
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs = 0.0_real64
        vecs(4:6, 1) = [5.0, 5.0, 5.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: constant vector")

        call assert_equal_array_real( &
            projections(:, 1), &
            [0.0_real64, 0.0_real64, 0.0_real64], &
            3, &
            1d-12, &
            "test_omics_field_RAP_projection_constant_vector: Expected zero vector" &
            )
    end subroutine test_constant_vector

    !> Test orthogonal vector
    subroutine test_orthogonal_vector()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = 0.0_real64
        vecs(1:3, 1) = [1.0, 0.0, -1.0]
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: orthogonal vector")
    end subroutine test_orthogonal_vector

    !> Test no axes selected (error)
    subroutine test_no_axes()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.false., .false., .false.]
        vecs_mask = [.true., .true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_true(ierr /= 0, "ierr should be nonzero for no axes selected")
    end subroutine test_no_axes

    !> Test no vectors selected (error)
    subroutine test_no_vectors()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.true., .true., .true.]
        vecs_mask = [.false., .false., .false.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_true(ierr /= 0, "ierr should be nonzero for no vectors selected")
    end subroutine test_no_vectors

    !> Test mixed selection
    subroutine test_mixed_selection()
        implicit none

        real(real64), dimension(6, 3) :: vecs
        logical :: axes_mask(3), vecs_mask(3)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [6, 3])
        axes_mask = [.true., .false., .true.]
        vecs_mask = [.true., .false., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: mixed selection")
    end subroutine test_mixed_selection

    subroutine omics_field_RAP_projection2(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections, ierr)
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_fields
            !! number of vectors per axis
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(int32), intent(in) :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_selected_axes, n_selected_fields), intent(out) :: projections
            !! projected vectors
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        print *, fields

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_fields, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_dimension_size(n_selected_fields, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(fields, size(fields, kind=int32), ierr, arg_pos=1_int32)
        if (count(fields_selection_mask, kind=int32) > n_selected_fields) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) > n_selected_axes) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)

        if (is_err(ierr)) return

        call omics_field_RAP_projection_helper(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections)
    end subroutine omics_field_RAP_projection2

    !> Test non-square vecs
    subroutine test_non_square_vecs()
        real(real64), dimension(4, 2) :: vecs
        logical :: axes_mask(4), vecs_mask(2)
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs
        real(real64), allocatable :: projections(:, :)

        vecs(:, 1) = [1.0, 2.0, 3.0, 4.0]
        vecs(:, 2) = [4.0, 3.0, 2.0, 1.0]
        axes_mask = [.true., .true., .true., .true.]
        vecs_mask = [.true., .true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 4, 2, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: non-square vecs")
    end subroutine test_non_square_vecs

    !> Test concrete example
    subroutine test_concrete_example()
        real(real64), dimension(6, 1) :: vecs
        logical :: axes_mask(3), vecs_mask(1)
        real(real64), allocatable :: projections(:, :)
        real(real64), dimension(3) :: expected_projection_vec
        integer(int32) :: ierr, n_selected_axes, n_selected_vecs

        vecs(:, 1) = [1.0, -3.0, 1.1, 3.0, 6.0, 2.2]
        expected_projection_vec = vecs(1:3, 1) - vecs(4:6, 1)
        expected_projection_vec = expected_projection_vec - sum(expected_projection_vec)/3

        axes_mask = [.true., .true., .true.]
        vecs_mask = [.true.]
        n_selected_axes = count(axes_mask)
        n_selected_vecs = count(vecs_mask)

        allocate (projections(n_selected_axes, n_selected_vecs))
        projections = 1

        call omics_field_RAP_projection2(vecs, 3, 1, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
        call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: concrete example")

        call assert_equal_array_real( &
            projections(:, 1), &
            expected_projection_vec, &
            3, &
            1d-12, &
            "test_omics_field_RAP_projection_constant_vector: Calculated projection doesn't match expected" &
            )
    end subroutine test_concrete_example

end module mod_test_rap_tools_omics_field_RAP_projection
