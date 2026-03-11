! filepath: test/run_tests.f90
program main
  use, intrinsic :: iso_fortran_env, only: int32
  use test_suite
  use mod_test_compute_edf, only: get_all_tests_compute_edf
  use mod_test_bst, only: get_all_tests_bst
  use mod_test_kd_tree, only: get_all_tests_kd_tree
  use mod_test_sorting, only: get_all_tests_sorting
  use mod_test_get_outliers, only: get_all_tests_get_outliers
  use mod_test_loess_smoothing, only: get_all_tests_loess_smoothing
  use mod_test_normalize_by_std_dev, only: get_all_tests_normalize_by_std_dev
  use mod_test_quantile_normalization, only: get_all_tests_quantile_normalization
  use mod_test_log2_transformation, only: get_all_tests_log2_transformation
  use mod_test_calc_tiss_avg, only: get_all_tests_tiss_avg
  use mod_test_calc_fchange, only: get_all_tests_calc_fchange
  use mod_test_euclidean_distance, only: get_all_tests_euclidean_distance
  use mod_test_rap_tools_omics_vector_RAP_projection, only: get_all_tests_rap_tools_omics_vector_RAP_projection
  use mod_test_rap_tools_omics_field_RAP_projection, only: get_all_tests_rap_tools_omics_field_RAP_projection
  use mod_test_clock_hand_angles, only: get_all_tests_clock_hand_angles
  use mod_test_relative_axis_contributions, only: get_all_tests_relative_axis_contributions
  use mod_test_tissue_versatility, only: get_all_tests_tissue_versatility
  use mod_test_tox_data, only: get_all_tests_tox_data
  use mod_test_normalization_pipeline, only: get_all_tests_normalization_pipeline
  use mod_test_shift_vectors, only: get_all_tests_shift_vectors
  use mod_test_gene_centroids, only: get_all_tests_gene_centroids
  use mod_test_conversions, only: get_all_tests_conversions
  use mod_test_arrays, only: get_all_tests_arrays
  use mod_test_paralog_analysis, only: get_all_tests_paralog_analysis
  use mod_test_trajectory_contribution_analysis, only: get_all_tests_trajectory_contribution_analysis
  use mod_test_trajectory_normalization, only: get_all_tests_trajectory_normalization
  use mod_test_normalization_unit_length, only: get_all_tests_normalization_unit_length
  use mod_test_clustering, only: get_all_tests_clustering
  use mod_test_data_integration, only: get_all_tests_data_integration
  implicit none

  integer :: nargs
  character(len=128) :: requested_suite
  character(len=512) :: test_list

  call initialize_suites()
  call add_suite("bst", get_all_tests_bst)
  call add_suite("compute_edf", get_all_tests_compute_edf)
  call add_suite("k-d-tree", get_all_tests_kd_tree)
  call add_suite("sorting", get_all_tests_sorting)
  call add_suite("get_outliers", get_all_tests_get_outliers)
  call add_suite("loess_smoothing", get_all_tests_loess_smoothing)
  call add_suite("normalization", get_all_tests_normalize_by_std_dev)
  call add_suite("quantile_normalization", get_all_tests_quantile_normalization)
  call add_suite("log2_transformation", get_all_tests_log2_transformation)
  call add_suite("calc_tiss_avg", get_all_tests_tiss_avg)
  call add_suite("calc_fchange", get_all_tests_calc_fchange)
  call add_suite("euclidean_distance", get_all_tests_euclidean_distance)
  call add_suite("rap_tools_omics_vector_RAP_projection", get_all_tests_rap_tools_omics_vector_RAP_projection)
  call add_suite("rap_tools_omics_field_RAP_projection", get_all_tests_rap_tools_omics_field_RAP_projection)
  call add_suite("clock_hand_angles", get_all_tests_clock_hand_angles)
  call add_suite("relative_axis_contributions", get_all_tests_relative_axis_contributions)
  call add_suite("tissue_versatility", get_all_tests_tissue_versatility)
  call add_suite("tox_data", get_all_tests_tox_data)
  call add_suite("normalization_pipeline", get_all_tests_normalization_pipeline)
  call add_suite("shift_vectors", get_all_tests_shift_vectors)
  call add_suite("arrays", get_all_tests_arrays)
  call add_suite("gene_centroids", get_all_tests_gene_centroids)
  call add_suite("conversions", get_all_tests_conversions)
  call add_suite("paralog_analysis", get_all_tests_paralog_analysis)
  call add_suite("trajectory_contribution_analysis", get_all_tests_trajectory_contribution_analysis)
  call add_suite("trajectory_normalization", get_all_tests_trajectory_normalization)
  call add_suite("normalization_unit_length", get_all_tests_normalization_unit_length)
  call add_suite("clustering", get_all_tests_clustering)
  call add_suite("data_integration", get_all_tests_data_integration)

  nargs = command_argument_count()

  select case (nargs)
    case (0)
      call run_all_suites()
    case (1)
      call get_command_argument(1, requested_suite)
      call run_suite_all(trim(requested_suite))
    case (2)
      call get_command_argument(1, requested_suite)
      call get_command_argument(2, test_list)
      call run_suite_named(trim(requested_suite), trim(test_list))
    case default
      print *, "Too many arguments."
      call print_usage()
      stop 1
  end select
  contains
    !> Print usage information.
    subroutine print_usage()
      integer :: i

      print *, "Usage:"
      print *, "  run_tests"
      print *, "  run_tests <suite>"
      print *, "  run_tests <suite> <test1,test2,...>"
      print *, ""
      print *, "Available test suites:"
      do i = 1, size(available_suites)
        print *, "  ", trim(available_suites(i)%name)
      end do
    end subroutine print_usage
 
    
end program main
    


  