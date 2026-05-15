#include "macros.h"

!> Module for managing data structures and geometric calculations 
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, &
						  ERR_ALLOC_FAIL, ERR_DIM_MISMATCH, ERR_INVALID_INPUT
    implicit none
    private
    public :: allocate_vectors, deallocate_vectors
    public :: allocate_labels,  deallocate_labels
    public :: calc_label_sphere_radius
    public :: build_label_kd_tree
	public :: calculate_labels_as_density

contains

    !> Allocates a vector array with error checking.
    subroutine allocate_vectors( vectors, n_dimensions, n_vectors, ierr )
        !| Allocatable array to be initialized 
        real( real64 ), allocatable, intent( inout ) :: vectors( :, : )
        !| Number of rows (dimensions) 
        integer( int32 ), intent( in ) :: n_dimensions
        !| Number of columns (vectors) 
        integer( int32 ), intent( in ) :: n_vectors
        !| Error code output 
        integer( int32 ), intent( out ) :: ierr

        call set_ok( ierr )

        call validate_dimension_size( n_dimensions, ierr )
        if ( is_err( ierr ) ) return

        call validate_dimension_size( n_vectors, ierr )
        if ( is_err( ierr ) ) return

        if ( allocated( vectors ) ) deallocate( vectors )

        M_ALLOCATE( vectors( n_dimensions, n_vectors ) )

    end subroutine allocate_vectors

    !> Safely deallocates a vector array if it is currently allocated. 
    subroutine deallocate_vectors( vectors )
        !| The 2D array to deallocate 
        real( real64 ), allocatable, intent( inout ) :: vectors( :, : )

        if ( allocated( vectors ) ) deallocate( vectors )

    end subroutine deallocate_vectors

    !> Allocates a label array for analysis.
    subroutine allocate_labels( label_values, label_dimension, n_vectors, ierr )
        !| Allocatable label array 
        real( real64 ), allocatable, intent( inout ) :: label_values( :, : )
        !| Dimensions of the labels 
        integer( int32 ), intent( in ) :: label_dimension
        !| Number of vectors 
        integer( int32 ), intent( in ) :: n_vectors
        !| Error code output 
        integer( int32 ), intent( out ) :: ierr

        call set_ok( ierr )

        call validate_dimension_size( label_dimension, ierr )
        if ( is_err( ierr ) ) return

        call validate_dimension_size( n_vectors, ierr )
        if ( is_err( ierr ) ) return

        if ( allocated( label_values ) ) deallocate( label_values )

        M_ALLOCATE( label_values( label_dimension, n_vectors ) )

    end subroutine allocate_labels

    !> Safely deallocates a label array. 
    subroutine deallocate_labels( label_values )
        !| The label array to deallocate 
        real( real64 ), allocatable, intent( inout ) :: label_values( :, : )

        if ( allocated( label_values ) ) deallocate( label_values )

    end subroutine deallocate_labels

    !> Calculates the radius of a label-sphere based on distance quantiles. 
    !! The radius is determined by finding the percentile of distances 
    !! from the mean vector to all other vectors in the set. 
    pure subroutine calc_label_sphere_radius( vectors, n_dimensions, n_vectors, &
                                              mean_vec, distances, perm, ierr, &
                                              radius, mean_to_other_vecs_dist_quant )

        use tox_gene_centroids, only: mean_vector
        use tox_euclidean_distance, only: euclidean_distance
        use f42_utils, only: sort_real_heapsort, calc_percentile

        !| Input data matrix (n_dimensions x n_vectors) 
        real( real64 ), intent( in ) :: vectors( :, : )
        !| Number of dimensions [cite: 139]
        integer( int32 ), intent( in ) :: n_dimensions
        !| Number of vectors [cite: 139]
        integer( int32 ), intent( in ) :: n_vectors
        !| Preallocated workspace for the mean vector 
        real( real64 ), intent( out ) :: mean_vec( : )
        !| Preallocated workspace for calculated distances 
        real( real64 ), intent( out ) :: distances( : )
        !| Preallocated workspace for permutation indices 
        integer( int32 ), intent( out ) :: perm( : )
        !| Error code output [cite: 139]
        integer( int32 ), intent( out ) :: ierr
        !| Optional quantile fraction (0.0 to 1.0), defaults to 0.15 
        real( real64 ), intent( in ), optional :: mean_to_other_vecs_dist_quant
        !| The resulting sphere radius 
        real( real64 ), intent( out ) :: radius

        real( real64 ) :: quantile_frac, percentile
        integer( int32 ) :: i

        call set_ok( ierr )

        call validate_dimension_size( n_dimensions, ierr )
        call validate_dimension_size( n_vectors, ierr )

        if ( is_err( ierr ) ) return

        if ( present( mean_to_other_vecs_dist_quant ) ) then
            quantile_frac = mean_to_other_vecs_dist_quant
        else
            quantile_frac = 0.15_real64
        end if

        call validate_in_range_real( quantile_frac, ierr, min = 0.0_real64, max = 1.0_real64 )
        if ( is_err( ierr ) ) return

        do i = 1, n_vectors
            perm( i ) = i
        end do

        call mean_vector( vectors, n_dimensions, n_vectors, perm, n_vectors, mean_vec, ierr )
        if ( is_err( ierr ) ) return

        do i = 1, n_vectors
            call euclidean_distance( mean_vec, vectors( :, i ), n_dimensions, distances( i ) )
        end do

        do i = 1, n_vectors
            perm( i ) = i
        end do

        call sort_real_heapsort( distances, perm )

        percentile = quantile_frac * 100.0_real64

        call calc_percentile( distances, perm, percentile, radius, ierr )

    end subroutine calc_label_sphere_radius

    !> Constructs KD-tree indices for the provided input vectors. 
    !! This routine acts as a wrapper for the core KD-tree indexing algorithm. 
    pure subroutine build_label_kd_tree( vectors, n_dimensions, n_vectors, dimension_order, kd_indices, &
                                         workspace, value_buffer, permutation, left_stack, right_stack, &
                                         recursion_stack, ierr )

        use f42_kd_tree, only: build_kd_index

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

        call validate_dimension_size( n_dimensions, ierr )
        call validate_dimension_size( n_vectors, ierr )

        if ( is_err( ierr ) ) return

        call build_kd_index( vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                             workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr )

    end subroutine build_label_kd_tree
	
	!> Density labels: count of vectors within radius r of each vector.
    !! Uses kd_indices ordering for traversal (no pruning yet).
    pure subroutine calculate_labels_as_density( vectors, n_dimensions, n_vectors, r, kd_indices, label_values, ierr )

        use tox_euclidean_distance, only: euclidean_distance

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