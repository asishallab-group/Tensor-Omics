! filepath: test/mod_test_field_vector_RAP_projection.f90
!> Unit test suite for RAP Projection routines.
module mod_test_rap_tools_omics_field_RAP_projection
   use asserts
   use tox_relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32
   use mod_test_suite, only: test_case
   implicit none
   public
   

contains

   !> Get array of all available tests.
   function get_all_tests_rap_tools_omics_field_RAP_projection() result(all_tests)
      type(test_case),allocatable :: all_tests(:)

      allocate(all_tests(10))

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

 
   
   
   !> Wrapper function for the actual call of `call_omics_field_RAP_projection`
   subroutine call_omics_field_RAP_projection(test_name, vecs, axes_mask, vecs_mask, result_projections)
      implicit none

      character(len=*), intent(in) :: test_name
      real(real64), dimension(:,:), intent(in) :: vecs
      logical, dimension(:), intent(in) :: axes_mask, vecs_mask
      real(real64), dimension(:,:), intent(out), optional :: result_projections

      integer(int32) :: n_axes, n_vecs, n_selected_axes, n_selected_vecs, ierr
      integer :: i_vec
      real(real64), allocatable :: projections(:,:)

      n_axes = size(axes_mask)
      n_vecs = size(vecs_mask)
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: " // trim(test_name))

      do i_vec = 1, n_selected_vecs
         call assert_equal_real(&
            sum(projections(:, i_vec)) / real(n_selected_axes, real64),&
            0.0_real64,&
            1d-12,&
            "test_omics_field_RAP_projection: Case '" // trim(test_name) // "' failed"&
         )
      end do

      if (present(result_projections)) then
         result_projections = projections
      end if
   end subroutine call_omics_field_RAP_projection

   !> Test all axes and vectors are selected
   subroutine test_all_selected()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs, i_vec
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: all selected")

      do i_vec = 1, n_selected_vecs
         call assert_equal_real(&
            sum(projections(:, i_vec)) / real(n_selected_axes, real64),&
            0.0_real64,&
            1d-12,&
            "test_omics_field_RAP_projection_all_selected: projection failed"&
         )
      end do
   end subroutine test_all_selected

   !> Test one axis and all vectors are selected
   subroutine test_one_axis_selected()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.false., .true., .false.]
      vecs_mask = [.true., .true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: one axis selected")
   end subroutine test_one_axis_selected

   !> Test all axes and one vector are selected
   subroutine test_one_vector_selected()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.false., .true., .false.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: one vector selected")
   end subroutine test_one_vector_selected

   !> Test constant vector
   subroutine test_constant_vector()
      implicit none

      real(real64), dimension(6,3) :: vecs
      real(real64), allocatable :: projections(:,:)
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs

      vecs = 0.0_real64
      vecs(4:6,1) = [5.0, 5.0, 5.0]
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: constant vector")

      call assert_equal_array_real(&
         projections(:, 1),&
         [0.0_real64, 0.0_real64, 0.0_real64],&
         3,&
         1d-12,&
         "test_omics_field_RAP_projection_constant_vector: Expected zero vector"&
      )
   end subroutine test_constant_vector

   !> Test orthogonal vector
   subroutine test_orthogonal_vector()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = 0.0_real64
      vecs(1:3,1) = [1.0, 0.0, -1.0]
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: orthogonal vector")
   end subroutine test_orthogonal_vector

   !> Test no axes selected (error)
   subroutine test_no_axes()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.false., .false., .false.]
      vecs_mask = [.true., .true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_true(ierr /= 0, "ierr should be nonzero for no axes selected")
   end subroutine test_no_axes

   !> Test no vectors selected (error)
   subroutine test_no_vectors()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.false., .false., .false.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_true(ierr /= 0, "ierr should be nonzero for no vectors selected")
   end subroutine test_no_vectors

   !> Test mixed selection
   subroutine test_mixed_selection()
      implicit none

      real(real64), dimension(6,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [6,3])
      axes_mask = [.true., .false., .true.]
      vecs_mask = [.true., .false., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 3, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: mixed selection")
   end subroutine test_mixed_selection

   !> Test non-square vecs
   subroutine test_non_square_vecs()
      real(real64), dimension(8,2) :: vecs
      logical :: axes_mask(4), vecs_mask(2)
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)

      vecs(1:4,1) = [1.0, 2.0, 3.0, 4.0]
      vecs(1:4,2) = [4.0, 3.0, 2.0, 1.0]
      axes_mask = [.true., .true., .true., .true.]
      vecs_mask = [.true., .true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 4, 2, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: non-square vecs")
   end subroutine test_non_square_vecs

   !> Test concrete example
   subroutine test_concrete_example()
      real(real64), dimension(6,1) :: vecs
      logical :: axes_mask(3), vecs_mask(1)
      real(real64), allocatable :: projections(:,:)
      real(real64), dimension(3) :: expected_projection_vec
      integer(int32) :: ierr, n_selected_axes, n_selected_vecs

      vecs(:,1) = [1.0, -3.0, 1.1, 3.0, 6.0, 2.2]
      expected_projection_vec = vecs(1:3,1) - vecs(4:6,1)
      expected_projection_vec = expected_projection_vec - sum(expected_projection_vec) / 3

      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true.]
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_field_RAP_projection(vecs, 3, 1, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections, ierr)
      call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: concrete example")

      call assert_equal_array_real(&
         projections(:, 1),&
         expected_projection_vec,&
         3,&
         1d-12,&
         "test_omics_field_RAP_projection_constant_vector: Calculated projection doesn't match expected"&
      )
   end subroutine test_concrete_example

end module mod_test_rap_tools_omics_field_RAP_projection