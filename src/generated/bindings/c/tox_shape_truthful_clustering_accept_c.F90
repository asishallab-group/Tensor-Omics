#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_accept(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_accept_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: accept_ensemble_expert_c
    public :: accept_ensemble_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept(module):accept_ensemble(subroutine)]]
    !| Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
    !| bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
    !| cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
    !| dimension itself changed, the two tangent bases don't share a common dimension to
    !| compare angles over at all, so this criterion is skipped (no SVD is computed) and
    !| treated as vacuously satisfied; criterion (2) is what actually judges whether that
    !| change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
    !| `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
    !| check could foresee.
    subroutine accept_ensemble_expert_c(&
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
        ) bind(C, name="accept_ensemble_expert_c")
        use tox_shape_truthful_clustering_accept, only: accept_ensemble

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: d_t
            !! Ensemble's intrinsic dimension at t
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_accept_kernel(module):tox_stc_accept_ensemble_svd_workspace]].
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_t
            !! Ensemble's tangent+normal basis at t, see `observable`
        real(c_double), intent(in), target :: G_t
            !! Ensemble's spectral gap at t
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_tp1
            !! Ensemble's tangent+normal basis at t+1
        real(c_double), intent(in), target :: G_tp1
            !! Ensemble's spectral gap at t+1
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: alpha_max
            !! Maximum tolerated principal angle (radians)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(min(d_t,d_tp1), min(d_t,d_tp1)), intent(out), target :: tmp_m
            !! Workspace: M = U_t(:,1:d)^T U_tp1(:,1:d)
        real(c_double), dimension(min(d_t,d_tp1)), intent(out), target :: tmp_s
            !! Workspace: singular values of M (= cos of the principal angles)
        real(c_double), dimension(lwork), intent(out), target :: tmp_work
            !! Workspace: LAPACK dgesvd scratch
        logical(c_bool), intent(out), target :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        logical :: is_accepted_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(d_t)
        M_CHECK_NON_NULL(G_t)
        M_CHECK_NON_NULL(d_tp1)
        M_CHECK_NON_NULL(G_tp1)
        M_CHECK_NON_NULL(alpha_max)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(lwork)
        M_CHECK_NON_NULL(is_accepted)
        M_CHECK_ARRAY_NON_NULL(U_t, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(U_tp1, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_m, (min(d_t,d_tp1)) * (min(d_t,d_tp1)))
        M_CHECK_ARRAY_NON_NULL(tmp_s, (min(d_t,d_tp1)))
        M_CHECK_ARRAY_NON_NULL(tmp_work, lwork)

        call accept_ensemble(&
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
            is_accepted = is_accepted_f,&
            ierr = ierr&
        )

        is_accepted = is_accepted_f
    end subroutine accept_ensemble_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_accept(module):accept_ensemble_alloc(subroutine)]]
    !| Three criteria, all must hold: (1) principal angles between the d-dimensional tangent
    !| bases, via `dgesvd` on M = U_t(:,1:d)^T U_tp1(:,1:d), whose singular values are
    !| cos(alpha_i) directly -- but only when d_t == d_tp1: when the estimated intrinsic
    !| dimension itself changed, the two tangent bases don't share a common dimension to
    !| compare angles over at all, so this criterion is skipped (no SVD is computed) and
    !| treated as vacuously satisfied; criterion (2) is what actually judges whether that
    !| change in d is acceptable. (2) |d_tp1 - d_t| <= d_max. (3) |log(G_tp1/G_t)| <= G_max.
    !| `ierr` is set only if the LAPACK SVD fails to converge -- not a condition any input
    !| check could foresee.
    subroutine accept_ensemble_c(&
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
        ) bind(C, name="accept_ensemble_c")
        use tox_shape_truthful_clustering_accept, only: accept_ensemble_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_t
            !! Ensemble's tangent+normal basis at t, see `observable`
        integer(c_int), intent(in), target :: d_t
            !! Ensemble's intrinsic dimension at t
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), target :: G_t
            !! Ensemble's spectral gap at t
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(in), target :: U_tp1
            !! Ensemble's tangent+normal basis at t+1
        integer(c_int), intent(in), target :: d_tp1
            !! Ensemble's intrinsic dimension at t+1
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), target :: G_tp1
            !! Ensemble's spectral gap at t+1
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: alpha_max
            !! Maximum tolerated principal angle (radians)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `2.0_real64 * atan(1.0_real64)`.
        integer(c_int), intent(in), target :: d_max
            !! Maximum tolerated change in intrinsic dimension
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(in), target :: G_max
            !! Maximum tolerated |log(G_tp1/G_t)|
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(out), target :: is_accepted
            !! .true. if all three acceptance criteria are satisfied
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        logical :: is_accepted_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(d_t)
        M_CHECK_NON_NULL(G_t)
        M_CHECK_NON_NULL(d_tp1)
        M_CHECK_NON_NULL(G_tp1)
        M_CHECK_NON_NULL(alpha_max)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(is_accepted)
        M_CHECK_ARRAY_NON_NULL(U_t, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(U_tp1, n_dimensions * n_dimensions)

        call accept_ensemble_alloc(&
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
            is_accepted = is_accepted_f,&
            ierr = ierr&
        )

        is_accepted = is_accepted_f
    end subroutine accept_ensemble_c

end module tox_shape_truthful_clustering_accept_c
#endif
