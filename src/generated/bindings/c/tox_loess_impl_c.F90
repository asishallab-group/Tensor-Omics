#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_loess_impl(module)]]
!| LOESS local polynomial regression smoothing, over the netlib `dloess`/`lowesd` family.
!|
!| Two fits: `loess_fit_plain` and `loess_fit_robust`, the latter reweighting against outliers over
!| a number of iterations. There is deliberately no combined entry point dispatching on a mode --
!| a caller chooses the routine, and supplies the weights and the evaluation points it wants.
!|
!| A sample too degenerate to fit (fewer distinct points than the degree needs) is answered
!| directly from the observations rather than handed to netlib, which cannot fit one and dies
!| inside its decomposition instead of reporting. Nothing about that is visible to a caller
!| beyond the result.
!|
!| `tox_loess_required_workspace` sizes the workspace for a caller that wants to own it; the
!| mode and iteration constants and `EPS_LOESS` are here for the same reason. The netlib
!| routines themselves are not re-documented beyond their calling convention.
module tox_loess_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: tox_loess_required_workspace_c

contains

    !> summary: C-wrapper for [[tox_loess_impl(module):tox_loess_required_workspace(subroutine)]]
    !| Computes the required sizes for integer and real workspace arrays.
    !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
    subroutine tox_loess_required_workspace_c(&
            n_dim,&
            max_neighborhood_size,&
            int_workspace_size,&
            real_workspace_size,&
            save_factorization,&
            ierr&
        ) bind(C, name="tox_loess_required_workspace_c")
        use tox_loess_impl, only: tox_loess_required_workspace

        integer(c_int), intent(in), target :: n_dim
            !! Dimensionality of the data
        integer(c_int), intent(in), target :: max_neighborhood_size
            !! Maximum neighborhood size
        integer(c_int), intent(out), target :: int_workspace_size
            !! Required size of the integer workspace array
        integer(c_int), intent(out), target :: real_workspace_size
            !! Required size of the real workspace array
        logical(c_bool), intent(in), target :: save_factorization
            !! Save matrix factorization flag
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dim)
        M_CHECK_NON_NULL(max_neighborhood_size)
        M_CHECK_NON_NULL(int_workspace_size)
        M_CHECK_NON_NULL(real_workspace_size)
        M_CHECK_NON_NULL(save_factorization)

        call tox_loess_required_workspace(&
            n_dim = n_dim,&
            max_neighborhood_size = max_neighborhood_size,&
            int_workspace_size = int_workspace_size,&
            real_workspace_size = real_workspace_size,&
            save_factorization = save_factorization&
        )
    end subroutine tox_loess_required_workspace_c

end module tox_loess_impl_c
#endif
