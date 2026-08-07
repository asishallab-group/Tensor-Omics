#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_observable_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_observable
    use tox_shape_truthful_clustering_observable_kernel, only: normal_error_kernel, observable_kernel, tangent_scales_kernel, tox_stc_observable_svd_workspace
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: normal_error
    public :: tangent_scales
    public :: observable
    public :: observable_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_kernel(module):normal_error_kernel]].
    !| No pass over the ensemble's member vectors is required; the sum is already implied by
    !| the singular value decomposition [[tox_shape_truthful_clustering_observable_kernel(module):observable_kernel]]
    !| computes.
    subroutine normal_error(&
            d,&
            eigenvalues,&
            n_dimensions,&
            normal_error_value,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: eigenvalues
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! The minimum valid value is `0.0_real64`.
        real(real64), intent(out) :: normal_error_value
            !! Mean squared residual off the d-dimensional tangent subspace
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(d, ierr, arg_pos=1_int32, min=0_int32, max=n_dimensions)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(eigenvalues, n_dimensions, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        call normal_error_kernel(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            normal_error_value = normal_error_value&
        )
    end subroutine normal_error

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_kernel(module):tangent_scales_kernel]].
    subroutine tangent_scales(&
            d,&
            eigenvalues,&
            n_dimensions,&
            tangent_scales_value,&
            ierr&
        )
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        real(real64), dimension(n_dimensions), intent(in) :: eigenvalues
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(d), intent(out) :: tangent_scales_value
            !! Extent along each of the d tangent directions
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(d, ierr, arg_pos=1_int32, min=0_int32, max=n_dimensions)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(eigenvalues, n_dimensions, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        call tangent_scales_kernel(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            tangent_scales_value = tangent_scales_value&
        )
    end subroutine tangent_scales

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_kernel(module):observable_kernel]].
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    subroutine observable(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            member_selection_mask,&
            n_selected_member,&
            lwork,&
            iwork_size,&
            tmp_y,&
            tmp_s,&
            tmp_u_econ,&
            tmp_vt_econ,&
            tmp_work,&
            tmp_iwork,&
            U,&
            eigenvalues,&
            mu,&
            d,&
            G,&
            normal_error_value,&
            tangent_scales_value,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in) :: n_selected_member
            !! Number of selected members (count of .TRUE. in member_selection_mask)
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_observable_kernel(module):tox_stc_observable_svd_workspace]].
        integer(int32), intent(in) :: iwork_size
            !! Size of tmp_iwork
            !! It is *VERY IMPORTANT* to compute this argument from the `iwork_size` output produced by [[tox_shape_truthful_clustering_observable_kernel(module):tox_stc_observable_svd_workspace]].
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        logical, dimension(n_vectors), intent(in) :: member_selection_mask
            !! Ensemble membership over the full dataset
        real(real64), dimension(n_dimensions, n_selected_member), intent(out) :: tmp_y
            !! Workspace: centered member matrix
        real(real64), dimension(min(n_dimensions,n_selected_member)), intent(out) :: tmp_s
            !! Workspace: singular values
        real(real64), dimension(n_dimensions, min(n_dimensions,n_selected_member)), intent(out) :: tmp_u_econ
            !! Workspace: economy-mode left singular vectors
        real(real64), dimension(min(n_dimensions,n_selected_member), n_selected_member), intent(out) :: tmp_vt_econ
            !! Workspace: economy-mode right singular vectors, transposed (unused beyond the SVD call)
        real(real64), dimension(lwork), intent(out) :: tmp_work
            !! Workspace: LAPACK dgesdd scratch
        integer(int32), dimension(iwork_size), intent(out) :: tmp_iwork
            !! Workspace: LAPACK dgesdd integer scratch
        real(real64), dimension(n_dimensions, n_dimensions), intent(out) :: U
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64), dimension(n_dimensions), intent(out) :: eigenvalues
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64), dimension(n_dimensions), intent(out) :: mu
            !! Ensemble center
        integer(int32), intent(out) :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64), intent(out) :: G
            !! Spectral gap at d
        real(real64), intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(real64), dimension(n_dimensions), intent(out) :: tangent_scales_value
            !! Extent along each tangent direction, zero-padded beyond d
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_member, ierr, arg_pos=5_int32, min=2_int32, max=n_vectors)
        call validate_dimension_size(lwork, ierr, arg_pos=6_int32)
        call validate_dimension_size(iwork_size, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (count(member_selection_mask, kind=int32) /= n_selected_member) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call observable_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            member_selection_mask = member_selection_mask,&
            n_selected_member = n_selected_member,&
            lwork = lwork,&
            iwork_size = iwork_size,&
            tmp_y = tmp_y,&
            tmp_s = tmp_s,&
            tmp_u_econ = tmp_u_econ,&
            tmp_vt_econ = tmp_vt_econ,&
            tmp_work = tmp_work,&
            tmp_iwork = tmp_iwork,&
            U = U,&
            eigenvalues = eigenvalues,&
            mu = mu,&
            d = d,&
            G = G,&
            normal_error_value = normal_error_value,&
            tangent_scales_value = tangent_scales_value,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine observable

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_observable_kernel(module):observable_kernel]].
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    subroutine observable_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            member_selection_mask,&
            n_selected_member,&
            U,&
            eigenvalues,&
            mu,&
            d,&
            G,&
            normal_error_value,&
            tangent_scales_value,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        logical, dimension(n_vectors), intent(in) :: member_selection_mask
            !! Ensemble membership over the full dataset
        integer(int32), intent(in) :: n_selected_member
            !! Number of selected members (count of .TRUE. in member_selection_mask)
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_dimensions, n_dimensions), intent(out) :: U
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64), dimension(n_dimensions), intent(out) :: eigenvalues
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64), dimension(n_dimensions), intent(out) :: mu
            !! Ensemble center
        integer(int32), intent(out) :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64), intent(out) :: G
            !! Spectral gap at d
        real(real64), intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(real64), dimension(n_dimensions), intent(out) :: tangent_scales_value
            !! Extent along each tangent direction, zero-padded beyond d
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success
        integer(int32) :: lwork
        integer(int32) :: iwork_size
        real(real64), dimension(:, :), allocatable :: tmp_y
        real(real64), dimension(:), allocatable :: tmp_s
        real(real64), dimension(:, :), allocatable :: tmp_u_econ
        real(real64), dimension(:, :), allocatable :: tmp_vt_econ
        real(real64), dimension(:), allocatable :: tmp_work
        integer(int32), dimension(:), allocatable :: tmp_iwork

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_member, ierr, arg_pos=5_int32, min=2_int32, max=n_vectors)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (count(member_selection_mask, kind=int32) /= n_selected_member) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call tox_stc_observable_svd_workspace(&
            n_dimensions = n_dimensions,&
            n_selected_member = n_selected_member,&
            lwork = lwork,&
            iwork_size = iwork_size&
        )
        M_ALLOCATE(tmp_y(n_dimensions, n_selected_member))
        M_ALLOCATE(tmp_s(min(n_dimensions,n_selected_member)))
        M_ALLOCATE(tmp_u_econ(n_dimensions, min(n_dimensions,n_selected_member)))
        M_ALLOCATE(tmp_vt_econ(min(n_dimensions,n_selected_member), n_selected_member))
        M_ALLOCATE(tmp_work(lwork))
        M_ALLOCATE(tmp_iwork(iwork_size))

        call observable_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            member_selection_mask = member_selection_mask,&
            n_selected_member = n_selected_member,&
            lwork = lwork,&
            iwork_size = iwork_size,&
            tmp_y = tmp_y,&
            tmp_s = tmp_s,&
            tmp_u_econ = tmp_u_econ,&
            tmp_vt_econ = tmp_vt_econ,&
            tmp_work = tmp_work,&
            tmp_iwork = tmp_iwork,&
            U = U,&
            eigenvalues = eigenvalues,&
            mu = mu,&
            d = d,&
            G = G,&
            normal_error_value = normal_error_value,&
            tangent_scales_value = tangent_scales_value,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine observable_alloc

end module tox_shape_truthful_clustering_observable
