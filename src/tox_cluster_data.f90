module tox_cluster_data
    
	use, intrinsic :: iso_fortran_env, only: real64
	implicit none
	
	type :: cluster_data
	    real( real64 ), allocatable ::  vectors( :, : )
		real( real64 ), allocatable ::  labels( :, : )
	end type cluster_data
	
end module tox_cluster_data