!> Unit test suite for normalization_pipeline routine.
module mod_test_normalization_pipeline
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains


  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(4)
    all_tests(1) = test_case("test_pipeline_basic", test_pipeline_basic)
    all_tests(2) = test_case("test_pipeline_edge_cases", test_pipeline_edge_cases)
    all_tests(3) = test_case("test_pipeline_vs_manual", test_pipeline_vs_manual)
    all_tests(4) = test_case("test_pipeline_empty_matrix", test_pipeline_empty_matrix)
  end function get_all_tests

  subroutine run_all_tests_normalization_pipeline()
    type(test_case) :: all_tests(4)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All normalization_pipeline tests passed successfully."
  end subroutine run_all_tests_normalization_pipeline

    !> Run specific normalization_pipeline tests by name.
  subroutine run_named_tests_normalization_pipeline(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(4)
    integer(int32) :: i, j
    logical :: found

    all_tests = get_all_tests()

    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_normalization_pipeline

  !> Basic test: pipeline with small matrix, no fold change
  subroutine test_pipeline_basic()
    integer(int32), parameter :: n_genes = 2, n_tissues = 2, n_grps = 2, max_stack = 10
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log
    real(real64), dimension(n_genes) :: temp_col, rank_means
    integer(int32), dimension(n_genes) :: perm
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    integer(int32) :: ierr

    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0] ! 2x2 matrix, col-major
    group_s = [1,2]; group_c = [1,1]

    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline returned error")
    call assert_no_nan_real(buf_log, n_genes * n_grps, "test_pipeline_basic: NaN in output")
    call assert_equal_int(size(buf_log), n_genes * n_grps, "test_pipeline_basic: output size incorrect")
    call assert_true(all(buf_log >= 0.0d0), "test_pipeline_basic: output should be non-negative")
  end subroutine test_pipeline_basic

  !> Edge case test: zero input, check output shape and values
  subroutine test_pipeline_edge_cases()
    integer(int32), parameter :: n_genes = 2, n_tissues = 2, n_grps = 2, max_stack = 10
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log
    real(real64), dimension(n_genes) :: temp_col, rank_means
    integer(int32), dimension(n_genes) :: perm
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    integer(int32) :: ierr

    input_matrix = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    group_s = [1,2]; group_c = [1,1]

    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline returned error")
    call assert_no_nan_real(buf_log, n_genes * n_grps, "test_pipeline_edge_cases: NaN in output")
    call assert_true(all(buf_log == 0.0d0), "test_pipeline_edge_cases: output should be all zeros for zero input")
  end subroutine test_pipeline_edge_cases

  !> Test: pipeline output matches manual stepwise normalization (no fold change)
  subroutine test_pipeline_vs_manual()
    integer(int32), parameter :: n_genes = 2, n_tissues = 2, n_grps = 2, max_stack = 10
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log, manual_out
    real(real64), dimension(n_genes) :: temp_col, rank_means
    integer(int32), dimension(n_genes) :: perm
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    integer(int32) :: ierr

    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
    group_s = [1,2]; group_c = [1,1]

    ! Manual stepwise normalization
    call normalize_by_std_dev(n_genes, n_tissues, input_matrix, buf_stddev, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev returned error")
    call quantile_normalization(n_genes, n_tissues, buf_stddev, buf_quant, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization returned error")
    call calc_tiss_avg(n_genes, n_grps, group_s, group_c, buf_quant, buf_avg, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg returned error")
    call log2_transformation(n_genes, n_tissues, buf_avg, manual_out, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation returned error")

    ! Pipeline normalization
    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline returned error")
    call assert_true(all(abs(buf_log - manual_out) < 1.0d-12), "test_pipeline_vs_manual: pipeline and manual outputs differ")
  end subroutine test_pipeline_vs_manual

  !> Test with empty input matrix.
  subroutine test_pipeline_empty_matrix()
    integer(int32) :: n_genes, n_tissues, n_grps, max_stack, ierr
    real(real64), dimension(0) :: input_matrix, buf_stddev, buf_quant, buf_avg, buf_log
    real(real64), dimension(0) :: temp_col, rank_means
    integer(int32), dimension(0) :: perm, stack_left, stack_right, group_s, group_c
    n_genes = 0; n_tissues = 0; n_grps = 0; max_stack = 0
    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, temp_col, rank_means, perm, stack_left, stack_right, max_stack, group_s, group_c, n_grps, ierr)
    call assert_equal_int(ierr, 202, "normalization_pipeline should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_pipeline_empty_matrix

end module mod_test_normalization_pipeline
