!> What the children of the `tox_clustering` suite share: the values an output holds before a
!| call, the three non-finite values every procedure has to refuse, and the check that a k-means
!| run reported as successful handed back an assignment its documentation allows.
module mod_test_tox_clustering_fixtures
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use asserts, only: assert_in_range_int, assert_equal_int, operator(//)
    implicit none
    private
    public :: UNWRITTEN_INT, UNWRITTEN_REAL, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    public :: check_documented_assignment

    !> What an integer output (a label, a count, a merge, a cluster size) holds before a call, so
    !| that an element the procedure never writes shows. No expected value in this suite equals it.
    integer(int32), parameter :: UNWRITTEN_INT = -777_int32
    !> What a real output (a height) holds before a call; no expected value equals it either.
    real(real64), parameter :: UNWRITTEN_REAL = -777.0_real64

    !> How many non-finite values `non_finite_value` hands out.
    integer(int32), parameter :: N_NON_FINITE = 3
    !> The names of the non-finite values, for the assertion messages.
    character(*), parameter :: NON_FINITE_NAMES(N_NON_FINITE) = ["NaN ", "+Inf", "-Inf"]

contains

    !> The `i_value`-th non-finite value: NaN, +Inf, -Inf.
    function non_finite_value(i_value) result(value)
        integer(int32), intent(in) :: i_value
        real(real64) :: value

        select case (i_value)
        case (1)
            value = ieee_value(1.0_real64, ieee_quiet_nan)
        case (2)
            value = ieee_value(1.0_real64, ieee_positive_inf)
        case default
            value = ieee_value(1.0_real64, ieee_negative_inf)
        end select
    end function non_finite_value

    !> What the k-means documentation promises of a run it reports as successful, whatever the
    !| centroids: every label lies in `[1, n_clusters]` ("`1<=label<=n_clusters=k`"), and
    !| `label_counts(k)` is the number of points with label `k`, so the counts add up to `n_points`.
    subroutine check_documented_assignment(n_points, n_clusters, labels, label_counts, msg)
        integer(int32), intent(in) :: n_points, n_clusters
        integer(int32), intent(in) :: labels(n_points), label_counts(n_clusters)
        character(*), intent(in) :: msg

        integer(int32) :: i_point, i_cluster

        do i_point = 1, n_points
            call assert_in_range_int(labels(i_point), 1_int32, n_clusters, msg//": label of point "//i_point)
        end do
        do i_cluster = 1, n_clusters
            call assert_equal_int(label_counts(i_cluster), count(labels == i_cluster, kind=int32), &
                                  msg//": label_counts("//i_cluster//") counts the points labelled "//i_cluster)
        end do
    end subroutine check_documented_assignment
end module mod_test_tox_clustering_fixtures
