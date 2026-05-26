#include "macros.h"

!> Module for managing data structures and geometric calculations 
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, &
						  ERR_ALLOC_FAIL, ERR_DIM_MISMATCH, ERR_INVALID_INPUT
	use tox_gene_centroids, only: mean_vector
	use tox_euclidean_distance, only: euclidean_distance
    use f42_utils, only: sort_real_heapsort, calc_percentile
	use f42_kd_tree, only: build_kd_index
    implicit none
    private
    public :: vectors_alloc
    public :: label_1d_alloc, labels_2d_alloc
    public :: calculate_label_sphere_radius
    public :: compute_label_kd_tree
	public :: calculate_labels_as_density

contains
	
    !> Convenience Allocating Wrapper for calculating label sphere radius.
	subroutine calculate_label_sphere_radius_alloc( vectors, n_dimensions, n_vectors, &
													radius, mean_to_other_vecs_dist_quant, ierr )

		integer( int32 ), intent( in ) :: n_dimensions
		!! Number of dimensions
		integer( int32 ), intent( in ) :: n_vectors
		!! Number of vectors
		real( real64 ), intent( in ) :: vectors( n_dimensions, n_vectors )
		!! Input data matrix (n_dimensions x n_vectors)
		real( real64 ), intent( out ) :: radius
		!! The resulting sphere radius
		real( real64 ), intent( in ), optional :: mean_to_other_vecs_dist_quant
		!! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
		integer( int32 ), intent( out ) :: ierr
		!! Error code output (always the final argument).

		! Local work arrays utilizing mandated tmp_ prefix
		real( real64 ), allocatable :: tmp_mean_vec(:)
		real( real64 ), allocatable :: tmp_distances(:)
		integer( int32 ), allocatable :: tmp_perm(:)

		call set_ok( ierr )

		! Perform structural allocations using the project safety macro
		M_ALLOCATE( tmp_mean_vec( n_dimensions ) )
		M_ALLOCATE( tmp_distances( n_vectors ) )
		M_ALLOCATE( tmp_perm( n_vectors ) )

		! Forward execution down to Layer 2 
		call calculate_label_sphere_radius( vectors, n_dimensions, n_vectors, &
											tmp_mean_vec, tmp_distances, tmp_perm, &
											radius, mean_to_other_vecs_dist_quant, ierr )

	end subroutine calculate_label_sphere_radius_alloc
		
	!> Validated Entry Point for calculating label sphere radius.
	subroutine calculate_label_sphere_radius( vectors, n_dimensions, n_vectors, &
											  mean_vec, distances, perm, radius, &
											  mean_to_other_vecs_dist_quant, ierr )

		integer( int32 ), intent( in ) :: n_dimensions
		!! Number of dimensions
		integer( int32 ), intent( in ) :: n_vectors
		!! Number of vectors
		real( real64 ), intent( in ) :: vectors( n_dimensions, n_vectors )
		!! Input data matrix (n_dimensions x n_vectors)
		real( real64 ), intent( out ) :: mean_vec( n_dimensions )
		!! Preallocated workspace for the mean vector
		real( real64 ), intent( out ) :: distances( n_vectors )
		!! Preallocated workspace for calculated distances
		integer( int32 ), intent( out ) :: perm( n_vectors )
		!! Preallocated workspace for permutation indices
		real( real64 ), intent( out ) :: radius
		!! The resulting sphere radius
		real( real64 ), intent( in ), optional :: mean_to_other_vecs_dist_quant
		!! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
		integer( int32 ), intent( out ) :: ierr
		!! Error code output (always the final argument).

		real( real64 ) :: actual_quant

		call set_ok( ierr )

		! Validate all input dimensions first
		call validate_dimension_size( n_dimensions, ierr )
		call validate_dimension_size( n_vectors, ierr )
		if ( is_err( ierr ) ) return

		! Assign default value using mandated preprocessor macro
		M_DEFAULT_VAL( mean_to_other_vecs_dist_quant, actual_quant, 0.15_real64 )
	
		! Validate the calculated/assigned range
		call validate_in_range_real( actual_quant, ierr, min = 0.0_real64, max = 1.0_real64 )
		if ( is_err( ierr ) ) return
	
		! Delegate directly to pure core implementation layer
		call calculate_label_sphere_radius_helper( vectors, n_dimensions, n_vectors, &
												   mean_vec, distances, perm, actual_quant, &
												   radius, ierr )
	
	end subroutine calculate_label_sphere_radius
	
	!> Core Implementation for calculating label sphere radius.
	pure subroutine calculate_label_sphere_radius_helper( vectors, n_dimensions, n_vectors, &
														  mean_vec, distances, perm, &
														  quantile_frac, radius, ierr )

		integer( int32 ), intent( in ) :: n_dimensions
		!! Number of dimensions
		integer( int32 ), intent( in ) :: n_vectors
		!! Number of vectors
		real( real64 ), intent( in ) :: vectors( n_dimensions, n_vectors )
		!! Input data matrix (n_dimensions x n_vectors)
		real( real64 ), intent( out ) :: mean_vec( n_dimensions )
		!! Preallocated workspace for the mean vector
		real( real64 ), intent( out ) :: distances( n_vectors )
		!! Preallocated workspace for calculated distances
		integer( int32 ), intent( out ) :: perm( n_vectors )
		!! Preallocated workspace for permutation indices
		real( real64 ), intent( in ) :: quantile_frac
		!! Quantile fraction verified (0.0 to 1.0)
		real( real64 ), intent( out ) :: radius
		!! The resulting sphere radius
		integer( int32 ), intent( out ) :: ierr
		!! Error code output (always the final argument).

		real( real64 ) :: percentile
		integer( int32 ) :: i_vec

		call set_ok( ierr )

		! Element-wise initialization executed in parallel
		do concurrent ( i_vec = 1:n_vectors ) shared( perm )
			perm( i_vec ) = i_vec
		end do

		call mean_vector( vectors, n_dimensions, n_vectors, perm, n_vectors, mean_vec, ierr )
		if ( is_err( ierr ) ) return

		! Multi-threaded distance execution via explicit locality specification
		do concurrent ( i_vec = 1:n_vectors ) shared( vectors, mean_vec, n_dimensions, distances )
			call euclidean_distance( mean_vec, vectors( :, i_vec ), n_dimensions, distances( i_vec ) )
		end do

		call sort_real_heapsort( distances, perm )

		percentile = quantile_frac * 100.0_real64

		call calc_percentile( distances, perm, percentile, radius, ierr )

	end subroutine calculate_label_sphere_radius_helper
		
    !> Constructs KD-tree indices for the provided input vectors. 
    !! This routine acts as a wrapper for the core KD-tree indexing algorithm. 
    pure subroutine compute_label_kd_tree( vectors, n_dimensions, n_vectors, dimension_order, kd_indices, &
                                         workspace, value_buffer, permutation, left_stack, right_stack, &
                                         recursion_stack, ierr )

        !| Input data matrix 
        real( real64 ), intent( in ) :: vectors( :, : )
        !| Number of dimensions 
        integer( int32 ), intent( in ) :: n_dimensions
        !| Number of vectors 
        integer( int32 ), intent( in ) :: n_vectors
        !| Order of dimensions for tree splitting 
        integer( int32 ), intent( in ) :: dimension_order( : )
        !| Output array for KD-tree indices 
        integer( int32 ), intent( out ) :: kd_indices( : )
        !| Integer workspace array 
        integer( int32 ), intent( out ) :: workspace( : )
        !| Real workspace buffer 
        real( real64 ), intent( out ) :: value_buffer( : )
        !| Permutation workspace 
        integer( int32 ), intent( out ) :: permutation( : )
        !| Workspace for the left branch stack 
        integer( int32 ), intent( out ) :: left_stack( : )
        !| Workspace for the right branch stack 
        integer( int32 ), intent( out ) :: right_stack( : )
        !| 2D array for recursion stack management 
        integer( int32 ), intent( out ) :: recursion_stack( :, : )
        !| Error code output 
        integer( int32 ), intent( out ) :: ierr

        call set_ok( ierr )
		
		! Validate dimensions and return immediately if an error occurs
        call validate_dimension_size( n_dimensions, ierr )
        if ( is_err( ierr ) ) return

        ! Validate vectors and return immediately if an error occurs
        call validate_dimension_size( n_vectors, ierr )
        if ( is_err( ierr ) ) return
		
        call build_kd_index( vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                             workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr )

    end subroutine compute_label_kd_tree
	
	!> Density labels: count of vectors within radius r of each vector.
    !! Uses kd_indices ordering for traversal (no pruning yet).
    pure subroutine calculate_labels_as_density( vectors, n_dimensions, n_vectors, r, kd_indices, label_values, ierr )

        !| Input vectors (n_dimensions x n_vectors), columns are vectors
        real( real64 ), intent( in ) :: vectors( :, : )
        !| Vector dimension ( rows )
        integer( int32 ), intent( in ) :: n_dimensions
        !| Number of vectors ( columns )
        integer( int32 ), intent( in ) :: n_vectors
        !| Label-sphere radius
        real( real64 ), intent( in ) :: r
        !| KD-tree index (order of vectors)
        integer( int32 ), intent( in ) :: kd_indices( : )
        !| Output labels (label_dimension x n_vectors), uses row 1
        real( real64 ), intent( out ) :: label_values( :, : )
        !| Error checker
        integer( int32 ), intent( out ) :: ierr

        integer( int32 ) :: i, k, idx, count
        real( real64 ) :: dist

        call set_ok( ierr )

        call validate_dimension_size( n_dimensions, ierr )
        if ( is_err( ierr ) ) return
        call validate_dimension_size( n_vectors, ierr )
        if ( is_err( ierr ) ) return

        call validate_in_range_real( r, ierr, min = 0.0_real64, max = huge( 1.0_real64 ) )
        if ( is_err( ierr ) ) return

        if ( size( vectors, 1 ) /= n_dimensions .or. size( vectors, 2 ) /= n_vectors ) then
            call set_err( ierr, ERR_DIM_MISMATCH )
            return
        end if

        if ( size( kd_indices ) < n_vectors ) then
            call set_err( ierr, ERR_DIM_MISMATCH )
            return
        end if

        if ( size( label_values, 2 ) /= n_vectors ) then
            call set_err( ierr, ERR_DIM_MISMATCH )
            return
        end if

        if ( size( label_values, 1 ) < 1 ) then
            call set_err( ierr, ERR_DIM_MISMATCH )
            return
        end if

        do i = 1, n_vectors
            label_values(:, i) = 0.0_real64
            count = 0

            do k = 1, n_vectors
                idx = kd_indices( k )
                if ( idx < 1 .or. idx > n_vectors ) then
                    call set_err( ierr, ERR_INVALID_INPUT )
                    return
                end if

                call euclidean_distance( vectors(:, i), vectors(:, idx), n_dimensions, dist )
                if ( dist <= r ) count = count + 1
            end do

            label_values(1, i) = real( count, real64 )
        end do

    end subroutine calculate_labels_as_density

end module tox_shatter_cluster_data
