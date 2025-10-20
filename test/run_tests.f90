program main
  use mod_test_bst
  use mod_test_kd_tree
  use mod_test_sorting
  use mod_test_get_outliers
  use mod_test_loess_smoothing
  use mod_test_normalize_by_std_dev
  use mod_test_quantile_normalization
  use mod_test_log2_transformation
  use mod_test_calc_tiss_avg
  use mod_test_calc_fchange
  use mod_test_euclidean_distance
  use mod_test_rap_tools_omics_vector_RAP_projection
  use mod_test_rap_tools_omics_field_RAP_projection
  use mod_test_clock_hand_angles
  use mod_test_relative_axis_contributions
  use mod_test_tissue_versatility
  use mod_test_normalization_pipeline
  use mod_test_shift_vectors
  use mod_test_gene_centroids
  use mod_test_tox_conversions
  use mod_test_arrays


  implicit none

  ! Type for suite registry
  type :: suite_entry
    character(len=64) :: name
    procedure(run_all_interface), pointer, nopass :: run_all => null()
    procedure(run_named_interface), pointer, nopass :: run_named => null()
  end type suite_entry

  ! Abstract interfaces
  abstract interface
    subroutine run_all_interface()
    end subroutine run_all_interface
    
    subroutine run_named_interface(test_names)
      character(len=*), intent(in) :: test_names(:)
    end subroutine run_named_interface
  end interface

  ! Registry of all available suites
  type(suite_entry), allocatable :: available_suites(:)

  integer :: nargs
  character(len=64) :: requested_suite, test_list

  ! Initialize the suite registry
  call initialize_suites()

  nargs = command_argument_count()
  
  if (nargs == 0) then
    ! Run all tests from all suites
    call run_all_suites()
    
  else if (nargs == 1) then
    ! Run all tests in specified suite
    call get_command_argument(1, requested_suite)
    call run_suite_all(trim(requested_suite))
    
  else if (nargs == 2) then
    ! Run specific tests in suite
    call get_command_argument(1, requested_suite)
    call get_command_argument(2, test_list)
    call run_suite_named(trim(requested_suite), test_list)
    
  else
    print *, "Too many arguments"
    call print_usage()
    stop 1
  end if

contains

  !> Initialize the suite registry - ADD NEW SUITES HERE (no numbers!)
  subroutine initialize_suites()
    ! Start with empty registry
    allocate(available_suites(0))
    
    call add_suite("bst", run_all_tests_bst, run_named_tests_bst)
    call add_suite("k-d-tree", run_all_tests_kd_tree, run_named_tests_kd_tree)
    call add_suite("sorting", run_all_tests_sorting, run_named_tests_sorting)
    call add_suite("get_outliers",run_all_tests_get_outliers, run_named_tests_get_outliers)
    call add_suite("loess_smoothing",run_all_tests_loess_smoothing, run_named_tests_loess_smoothing)
    call add_suite("normalization", run_all_tests_normalize_by_std_dev, run_named_tests_normalize_by_std_dev)
    call add_suite("quantile_normalization", run_all_tests_quantile_normalization, run_named_tests_quantile_normalization)
    call add_suite("log2_transformation", run_all_tests_log2_transformation, run_named_tests_log2_transformation)
    call add_suite("calc_tiss_avg", run_all_tests_calc_tiss_avg, run_named_tests_calc_tiss_avg)
    call add_suite("calc_fchange", run_all_tests_calc_fchange, run_named_tests_calc_fchange)
    call add_suite("euclidean_distance", run_all_tests_euclidean_distance, run_named_tests_euclidean_distance)
    call add_suite("rap_tools_omics_vector_RAP_projection", run_all_tests_rap_tools_omics_vector_RAP_projection, run_named_tests_rap_tools_omics_vector_RAP_projection)
    call add_suite("rap_tools_omics_field_RAP_projection", run_all_tests_rap_tools_omics_field_RAP_projection, run_named_tests_rap_tools_omics_field_RAP_projection)
    call add_suite("clock_hand_angles", run_all_tests_clock_hand_angles, run_named_tests_clock_hand_angles)
    call add_suite("relative_axis_contributions", run_all_tests_relative_axis, run_named_tests_relative_axis)
    call add_suite("tissue_versatility", run_all_tests_tissue_versatility, run_named_tests_tissue_versatility)
    call add_suite("normalization_pipeline", run_all_tests_normalization_pipeline, run_named_tests_normalization_pipeline)
    call add_suite("shift_vectors", run_all_tests_shift_vectors, run_named_tests_shift_vectors)
    call add_suite("arrays", run_all_tests_array, run_named_tests_array)
    call add_suite("gene_centroids", run_all_tests_gene_centroids, run_named_tests_gene_centroids)
    call add_suite("tox_conversions", run_all_tests_tox_conversions, run_named_tests_tox_conversions)
    
  end subroutine initialize_suites
  

  !> Add a suite to the registry (grows automatically)
  subroutine add_suite(name, run_all_proc, run_named_proc)
    character(len=*), intent(in) :: name
    procedure(run_all_interface) :: run_all_proc
    procedure(run_named_interface) :: run_named_proc
    type(suite_entry), allocatable :: temp_suites(:)
    integer :: n
    
    n = size(available_suites)
    
    ! Create temporary array with one more slot
    allocate(temp_suites(n + 1))
    
    ! Copy existing suites
    if (n > 0) then
      temp_suites(1:n) = available_suites(1:n)
    end if
    
    ! Add new suite
    temp_suites(n + 1) = suite_entry(name, run_all_proc, run_named_proc)
    
    ! Replace the registry
    call move_alloc(temp_suites, available_suites)
  end subroutine add_suite

  !> Run all tests from all suites
  subroutine run_all_suites()
    integer :: i
    do i = 1, size(available_suites)
      print *, "Running suite: ", trim(available_suites(i)%name)
      call available_suites(i)%run_all()
    end do
  end subroutine run_all_suites

  !> Run all tests in a specific suite
  subroutine run_suite_all(requested_suite)
    character(len=*), intent(in) :: requested_suite
    integer :: i
    
    do i = 1, size(available_suites)
      if (trim(available_suites(i)%name) == requested_suite) then
        call available_suites(i)%run_all()
        return
      end if
    end do
    
    print *, "Unknown test suite: ", requested_suite
    call print_usage()
    stop 1
  end subroutine run_suite_all

  !> Run named tests in a specific suite
  subroutine run_suite_named(requested_suite, test_list)
    character(len=*), intent(in) :: requested_suite, test_list
    integer :: i
    
    do i = 1, size(available_suites)
      if (trim(available_suites(i)%name) == requested_suite) then
        call run_tests_from_list(test_list, available_suites(i)%run_named)
        return
      end if
    end do
    
    print *, "Unknown test suite: ", requested_suite
    call print_usage()
    stop 1
  end subroutine run_suite_named

  !> Run tests from comma-separated list using a specific runner
  subroutine run_tests_from_list(test_list, run_named_proc)
    character(len=*), intent(in) :: test_list
    procedure(run_named_interface) :: run_named_proc
    character(len=64) :: test_name
    character(len=64) :: single_test_array(1)
    integer :: start, end, pos
    
    start = 1
    
    do while (start <= len_trim(test_list))
      pos = index(test_list(start:), ',')
      if (pos > 0) then
        end = start + pos - 2
      else
        end = len_trim(test_list)
      end if
      
      test_name = trim(adjustl(test_list(start:end)))
      single_test_array(1) = test_name
      call run_named_proc(single_test_array)
      
      start = end + 2
      if (pos == 0) exit
    end do
  end subroutine run_tests_from_list

  !> Print usage information
  subroutine print_usage()
    integer :: i
    print *, "Usage:"
    print *, "  run_tests                                   # Run all tests"
    print *, "  run_tests <suite>                           # Run all tests in suite"
    print *, "  run_tests <suite> <test1,test2,...>         # Run specific tests"
    print *, ""
    print *, "Available test suites:"
    do i = 1, size(available_suites)
      print *, "  ", trim(available_suites(i)%name)
    end do
  end subroutine print_usage

end program main
