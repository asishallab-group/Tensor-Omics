#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_accept_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_accept
    use tox_shape_truthful_clustering_accept_kernel, only: accept_ensemble_kernel, tox_stc_accept_ensemble_svd_workspace
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math, only: above
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: accept_ensemble
    public :: accept_ensemble_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_accept_kernel(module):accept_ensemble_kernel]].
    !| Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
    !| bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
    !| cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
    !| dimension itself changed, the two tangent bases don't share a common dimension to
    !| compare angles over at all, so this criterion is skipped (no SVD is computed) and
    !| treated as vacuously satisfied; criterion (2) is what actually judges whether that
    !| change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
    !| `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
    !| check could foresee.
    subroutine accept_ensemble(&
            n_dimensions,&
            U_t,&
            d_t,&
            G_t,&
            U_tp1,&
            d_tp1,&
            G_tp1,&
            alpha_max,&
            d_max,&
            G_max,&
            lwork,&
            tmp_m,&
            tmp_s,&
            tmp_work,&
            is_accepted,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d_t
            !! Ensemble's intrinsic dimension at t
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_accept_kernel(module):tox_stc_accept_ensemble_svd_workspace]].
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_t
            !! Ensemble's tangent+normal basis at t, see `observable`
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at t
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_tp1
            !! Ensemble's tangent+normal basis at t+1
        real(real64), intent(in) :: G_tp1
            !! Ensemble's spectral gap at t+1
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(min(d_t,d_tp1), min(d_t,d_tp1)), intent(out) :: tmp_m
            !! Workspace: M = U_t(:,1:d)^T U_tp1(:,1:d)
        real(real64), dimension(min(d_t,d_tp1)), intent(out) :: tmp_s
            !! Workspace: singular values of M (= cos of the principal angles)
        real(real64), dimension(lwork), intent(out) :: tmp_work
            !! Workspace: LAPACK dgesvd scratch
        logical, intent(out) :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=1_int32)
        call validate_in_range_int(d_t, ierr, arg_pos=3_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_t, ierr, arg_pos=4_int32, min=above(0.0_real64))
        call validate_in_range_int(d_tp1, ierr, arg_pos=6_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_tp1, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_real(alpha_max, ierr, arg_pos=8_int32, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, arg_pos=9_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_dimension_size(lwork, ierr, arg_pos=11_int32)
        call validate_all_in_range_real(U_t, n_dimensions * n_dimensions, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(U_tp1, n_dimensions * n_dimensions, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call accept_ensemble_kernel(&
            n_dimensions = n_dimensions,&
            U_t = U_t,&
            d_t = d_t,&
            G_t = G_t,&
            U_tp1 = U_tp1,&
            d_tp1 = d_tp1,&
            G_tp1 = G_tp1,&
            alpha_max = alpha_max,&
            d_max = d_max,&
            G_max = G_max,&
            lwork = lwork,&
            tmp_m = tmp_m,&
            tmp_s = tmp_s,&
            tmp_work = tmp_work,&
            is_accepted = is_accepted,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine accept_ensemble

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_accept_kernel(module):accept_ensemble_kernel]].
    !| Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
    !| bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
    !| cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
    !| dimension itself changed, the two tangent bases don't share a common dimension to
    !| compare angles over at all, so this criterion is skipped (no SVD is computed) and
    !| treated as vacuously satisfied; criterion (2) is what actually judges whether that
    !| change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
    !| `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
    !| check could foresee.
    subroutine accept_ensemble_alloc(&
            n_dimensions,&
            U_t,&
            d_t,&
            G_t,&
            U_tp1,&
            d_tp1,&
            G_tp1,&
            alpha_max,&
            d_max,&
            G_max,&
            is_accepted,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_t
            !! Ensemble's tangent+normal basis at t, see `observable`
        integer(int32), intent(in) :: d_t
            !! Ensemble's intrinsic dimension at t
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), intent(in) :: G_t
            !! Ensemble's spectral gap at t
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(in) :: U_tp1
            !! Ensemble's tangent+normal basis at t+1
        integer(int32), intent(in) :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), intent(in) :: G_tp1
            !! Ensemble's spectral gap at t+1
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: alpha_max
            !! Maximum tolerated principal angle (radians)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(int32), intent(in) :: d_max
            !! Maximum tolerated change in intrinsic dimension
            !! The minimum valid value is `0_int32`.
        real(real64), intent(in) :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        logical, intent(out) :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success
        integer(int32) :: lwork
        real(real64), dimension(:, :), allocatable :: tmp_m
        real(real64), dimension(:), allocatable :: tmp_s
        real(real64), dimension(:), allocatable :: tmp_work

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=1_int32)
        call validate_in_range_int(d_t, ierr, arg_pos=3_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_t, ierr, arg_pos=4_int32, min=above(0.0_real64))
        call validate_in_range_int(d_tp1, ierr, arg_pos=6_int32, min=0_int32, max=n_dimensions)
        call validate_in_range_real(G_tp1, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_real(alpha_max, ierr, arg_pos=8_int32, min=0.0_real64, max=2.0_real64 * atan(1.0_real64))
        call validate_in_range_int(d_max, ierr, arg_pos=9_int32, min=0_int32)
        call validate_in_range_real(G_max, ierr, arg_pos=10_int32, min=0.0_real64)
        call validate_all_in_range_real(U_t, n_dimensions * n_dimensions, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(U_tp1, n_dimensions * n_dimensions, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call tox_stc_accept_ensemble_svd_workspace(&
            d_t = d_t,&
            d_tp1 = d_tp1,&
            lwork = lwork&
        )
        M_ALLOCATE(tmp_m(min(d_t,d_tp1), min(d_t,d_tp1)))
        M_ALLOCATE(tmp_s(min(d_t,d_tp1)))
        M_ALLOCATE(tmp_work(lwork))

        call accept_ensemble_kernel(&
            n_dimensions = n_dimensions,&
            U_t = U_t,&
            d_t = d_t,&
            G_t = G_t,&
            U_tp1 = U_tp1,&
            d_tp1 = d_tp1,&
            G_tp1 = G_tp1,&
            alpha_max = alpha_max,&
            d_max = d_max,&
            G_max = G_max,&
            lwork = lwork,&
            tmp_m = tmp_m,&
            tmp_s = tmp_s,&
            tmp_work = tmp_work,&
            is_accepted = is_accepted,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine accept_ensemble_alloc

end module tox_shape_truthful_clustering_accept
