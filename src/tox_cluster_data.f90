!> Data structure for shape-truthful clustering
!!
!! Stores input vectors and their corresponding label values.
!!
!! Example:
!!
!! ```fortran
!! type( cluster_data ) :: c_data
!! allocate( c_data%vectors( 3, 75 ) )    ! 3D space, 75 vectors
!! allocate( c_data%labels( 1, 75 ) )     ! scalar density labels
!! c_data%n_dims = 3
!! c_data%n_vectors = 100
!! c_data%n_labels_dim = 1
!! ```
!!
module tox_cluster_data
    
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
	
    !> Data structure for clustering input vectors and their labels
    !!
    !! Stores input vectors in column-major format
    !! and associated label values with 1-to-1 column alignment.
    type :: cluster_data
    
        ! Number of dimensions
        integer( int32 ) :: n_dims
        ! Number of vectors        
        integer( int32 ) :: n_vectors 
        ! Dimensionality of labels    
        integer( int32 ) :: n_labels_dim 
        
        ! Input vectors: shape (n_dims, n_vectors)
        real( real64 ), allocatable ::  vectors( :, : ) 
        ! Label values: shape (n_labels_dim, n_vectors)
        real( real64 ), allocatable ::  labels( :, : ) 
        
    end type cluster_data
    
    contains

        !> Compute the mean (centroid) of all input vectors
        !! For each dimension, calculates the average across all vectors
		!!
		!! **Memory Allocation:**
		!! - Caller must preallocate `gene_indices` with size >= n_vectors
		!! - Caller must preallocate `mean_vec` with size = n_dims
		!! @note All input vectors are automatically selected; use mean_vector() directly
        pure subroutine compute_mean_vector( vectors, n_dims, n_vectors, gene_indices, mean_vec, ierr )
        
            use tox_gene_centroids, only: mean_vector
            use tox_errors, only: set_ok, set_err_once, ERR_EMPTY_INPUT, ERR_INVALID_INPUT
            
			! Number of dimensions
            integer( int32 ), intent( in ) :: n_dims
			! Number of vectors
            integer( int32 ), intent( in ) :: n_vectors
			! Input vectors matrix( n_dims, n_vectors ) -> Each column is a vectors.
            real( real64 ), intent( in ) :: vectors( :, : )
			! Preallocated arry to store indices
			integer( int32 ), intent( out) :: gene_indices( : )
			!Output mean vector
            real( real64 ), intent( out ) :: mean_vec( : )
			!Error checker
            integer( int32 ), intent( out ) :: ierr
            
            integer( int32 ) :: i
            
            call set_ok( ierr )
          
			! Validate inputs
			if (n_dims <= 0 .or. n_vectors <= 0) then
			
				call set_err_once(ierr, ERR_EMPTY_INPUT)
				
				return
				
			end if
			
			! Validate preallocated array sizes
			if ( size( gene_indices ) < n_vectors ) then
			
				call set_err_once( ierr, ERR_INVALID_INPUT )
				
				return
				
			end if

			if ( size( mean_vec ) < n_dims ) then
			
				call set_err_once( ierr, ERR_INVALID_INPUT )
				
				return
				
			end if
			
			! Create indices array
            do i = 1, n_vectors
            
                gene_indices( i ) = i
                
            end do
            
			
            call mean_vector( vectors, n_dims, n_vectors, gene_indices( 1: n_vectors ), n_vectors, &
							  mean_vec( 1:n_dims ), ierr )
            
        end subroutine compute_mean_vector
        
        !> Calculate Euclidean distances from mean vector to all input vectors
        !! For each vector, computes: sqrt(sum((vector - mean)^2))
        pure subroutine calculate_distances_from_mean( vectors, mean_vec, n_dims, n_vectors, distances, ierr )
        
            use tox_euclidean_distance, only: euclidean_distance
            use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
            
			! Number of dimensions
            integer( int32 ), intent( in ) :: n_dims
			! Number of vectors
            integer( int32 ), intent( in ) :: n_vectors
			! Input vectors matrix( n_dims, n_vectors ) -> Each column is a vectors.
            real( real64 ), intent( in ) :: vectors( :, : )
			! Mean vector
            real( real64 ), intent( in ) :: mean_vec( : )
			! Output distance vector
            real( real64 ), intent( out ) :: distances( : )
			! Error checker
            integer( int32 ), intent( out ) :: ierr
            
            integer( int32 ) :: i_vec
           
            call set_ok( ierr )
            
			! Validate inputs
            if ( n_dims <= 0 .or. n_vectors <= 0 ) then
            
                call set_err_once( ierr, ERR_INVALID_INPUT )
                
                distances = 0.0_real64
                
                return
                
            end if
            ! Validate array sizes match dimensions
			if ( size( vectors, 1 ) /= n_dims .or. size( vectors, 2 ) /= n_vectors ) then
			
				call set_err_once( ierr, ERR_INVALID_INPUT )
				
				distances = 0.0_real64
				
				return
				
			end if
			
			if ( size( mean_vec ) /= n_dims .or. size( distances ) < n_vectors ) then
			
				call set_err_once( ierr, ERR_INVALID_INPUT )
				
				distances = 0.0_real64
				
				return
				
			end if
			
			! Parallel loop: each thread calculates distances independently
            do i_vec = 1, n_vectors
			
				call euclidean_distance( vectors( :, i_vec ), mean_vec, n_dims, distances( i_vec ) )
				
			end do
            
        end subroutine calculate_distances_from_mean
        
        !> Get quantile value from sorted array of distances
        !! Sorts distances and returns the value at the quantile position
		!!
		!! @Quantile must be in range [0.0, 1.0].
		!! @note Caller must preallocate perm(:) with size >= n_distances
        pure subroutine sort_get_quantile_value( distances, n_distances, quantile, perm, q_value, ierr ) 
        
            use f42_utils, only: sort_real_heapsort
            use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
            
			! Number of distances
            integer( int32 ), intent( in ) :: n_distances
			! Input distances array unsorted
            real( real64 ), intent( in ) :: distances( : )
			! Desired quantile value
            real( real64 ), intent( in ) :: quantile
			! Permutation working array
			integer( int32 ), intent( out ) :: perm( : )
			! Output quantile value
            real( real64 ), intent( out ) :: q_value
			! Error checker
            integer( int32 ), intent( out ) :: ierr
			
            integer( int32 ) :: index, i
            
            call set_ok( ierr )
            
			! Validate inputs
            if ( n_distances <= 0 ) then
            
                call set_err_once( ierr, ERR_INVALID_INPUT )
                
                q_value = 0.0_real64
                
                return
                
            end if
            
            if ( quantile < 0.0_real64 .or. quantile > 1.0_real64 ) then
            
                call set_err_once( ierr, ERR_INVALID_INPUT )
                
                q_value = 0.0_real64
                
                return
                
            end if
            
			! Initialize array
            do i = 1, n_distances
            
                perm( i ) = i
                
            end do
            
			! Sort array.
            call sort_real_heapsort( distances( 1:n_distances ), perm( 1:n_distances ) )
            
			! Calculate position
            index = ceiling( real( quantile * n_distances, real64 ) )
            
            if ( index < 1 ) index = 1
            if ( index > n_distances ) index = n_distances
            
			! Quantile position
            q_value = distances( perm( index ) )
            
        end subroutine sort_get_quantile_value
        
        !> Calculate label-sphere radius using mean vector and distance quantile
        !!
        !! Algorithm:
        !! 1. Compute mean of all input vectors
        !! 2. Calculate distance from mean to each vector
        !! 3. Return the specified quantile of these distances as radius
        !!
        !! All arrays must be preallocated by the caller (per F42 standards).
        !!
        !! @note
        !! Workspace arrays (mean_vec, distances) are modified by this subroutine.
        !! The caller is responsible for allocating these with correct dimensions.
        !! @endnote
        pure subroutine calculate_radius_vectors( vectors, n_dims, n_vectors, &
                                                  mean_to_other_vecs_dist_quant, &
												  gene_indices, perm, &
                                                  mean_vec, distances, radius, ierr )
            
            use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT, is_err
            
			! Number of dimensions
            integer( int32 ), intent( in ) :: n_dims
			! Number of vectors
            integer( int32 ), intent( in ) :: n_vectors
			! Input vectors
            real( real64 ), intent( in ) :: vectors( :, : )
			! Quantile parameter
            real( real64 ), intent( in ) :: mean_to_other_vecs_dist_quant
			! Working array for indices 
			integer( int32 ), intent( out ) :: gene_indices( : )		
			! Working array for permutation
			integer( int32 ), intent( out ) :: perm( : )		
			! Output mean vector
            real( real64 ), intent( out ) :: mean_vec( : )
			!Output distance
            real( real64 ), intent( out ) :: distances( : )
			! Sphere radius
            real( real64 ), intent( out ) :: radius
			! Error checker
            integer( int32 ), intent( out ) :: ierr
            
            ! Initialize error
            call set_ok( ierr )
            
            ! Validate input dimensions
            if ( n_dims <= 0 .or. n_vectors <= 0 ) then
            
                call set_err_once( ierr, ERR_INVALID_INPUT )
                
                radius = 0.0_real64
                
                return
                
            end if
            
            ! Validate quantile parameter
            if ( mean_to_other_vecs_dist_quant < 0.0_real64 .or. &
                 mean_to_other_vecs_dist_quant > 1.0_real64 ) then
            
                call set_err_once( ierr, ERR_INVALID_INPUT )
                
                radius = 0.0_real64
                
                return
                
            end if
            
            ! Step 1: Compute mean vector
            call compute_mean_vector( vectors, n_dims, n_vectors, gene_indices(1:n_vectors), &
									  mean_vec(1:n_dims), ierr )
            
            if ( is_err( ierr ) ) then
            
                radius = 0.0_real64
                
                return
                
            end if
            
            ! Step 2: Calculate distances from mean to all vectors
            call calculate_distances_from_mean( vectors, mean_vec( 1:n_dims ), n_dims, n_vectors, &
                                                distances( 1:n_vectors ), ierr )
            
            if ( is_err( ierr ) ) then
            
                radius = 0.0_real64
                
                return
                
            end if
            
            ! Step 3: Get the quantile value as radius
            call sort_get_quantile_value( distances( 1:n_vectors ), n_vectors, &
                                          mean_to_other_vecs_dist_quant, perm( 1:n_vectors ), &
										  radius, ierr )
										  
        end subroutine calculate_radius_vectors
        
        !> Calculate scalar density labels for each vector.
        !!
        !! Density = count of vectors within sphere radius r
        !!
        !! @note
        !! All arrays must be preallocated by the caller.
        !! @endnote
        pure subroutine calculate_labels_as_density( vectors, n_dims, n_vectors, &
                                                     mean_to_other_vecs_dist_quant, &
                                                     mean_vec, distances, radius, labels, ierr )
            
            use tox_errors, only: set_ok
            
            integer( int32 ), intent( in ) :: n_dims
            integer( int32 ), intent( in ) :: n_vectors
            real( real64 ), intent( in ) :: vectors( n_dims, n_vectors )
            real( real64 ), intent( in ) :: mean_to_other_vecs_dist_quant
            real( real64 ), intent( out ) :: mean_vec( n_dims )
            real( real64 ), intent( out ) :: distances( n_vectors )
            real( real64 ), intent( out ) :: radius
            real( real64 ), intent( out ) :: labels( :, : )
            integer( int32 ), intent( out ) :: ierr
            
            ! Initialize error
            call set_ok( ierr )
            
            ! TODO: Implementation
            ! 1. Call calculate_radius_vectors to get radius
            ! 2. For each vector, count neighbors within radius
            ! 3. Store count in labels(:, i)
            
        end subroutine calculate_labels_as_density
        
end module tox_cluster_data
