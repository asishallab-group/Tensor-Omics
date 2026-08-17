#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Observable
!|
!| `observable`: the tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble,
!| obtained from the economy-mode singular value decomposition (LAPACK `dgesdd`) of its
!| centered member vectors -- never an eigendecomposition of an explicitly formed
!| covariance matrix (see `misc/mod_STC.md`, "Numerical Linear Algebra"). `normal_error` and
!| `tangent_scales` are simple, dependency-free reductions over the eigenvalues `observable`
!| computes. See `misc/mod_STC.md`, SKG `observable`/`normal_error`/`tangent_scales`, for the
!| full algorithm definitions.
!|
!| Generated from [[tox_shape_truthful_clustering_observable_impl(module)]]; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_observable
    use f42_safeguard
    use tox_shape_truthful_clustering_observable_impl, only: ensemble_final_observable_impl, normal_error_impl, observable_impl, tangent_scales_impl
    use tox_shape_truthful_clustering_observable_impl, only: tox_stc_observable_svd_workspace
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_real
    use tox_errors, only: validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: normal_error
    public :: tangent_scales
    public :: observable
    public :: observable_expert
    public :: ensemble_final_observable

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_impl(module):normal_error_impl]].
    !| No pass over the ensemble's member vectors is required; the sum is already implied by
    !| the singular value decomposition [[tox_shape_truthful_clustering_observable_impl(module):observable_impl]]
    !| computes.
    pure subroutine normal_error(&
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

        call normal_error_impl(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            normal_error_value = normal_error_value&
        )
    end subroutine normal_error

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_impl(module):tangent_scales_impl]].
    pure subroutine tangent_scales(&
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

        call tangent_scales_impl(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            tangent_scales_value = tangent_scales_value&
        )
    end subroutine tangent_scales

    !> summary: Validates its inputs, prepares what [[tox_shape_truthful_clustering_observable_impl(module):observable_impl]] needs, then calls it. The entry point to reach for first; see [[tox_shape_truthful_clustering_observable(module):observable_expert]] to prepare it yourself.
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    pure subroutine observable(&
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
        logical(c_bool), dimension(n_vectors), intent(in) :: member_selection_mask
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

        call observable_impl(&
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

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_impl(module):observable_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_shape_truthful_clustering_observable(module):observable]] does both.
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    pure subroutine observable_expert(&
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
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_observable_impl(module):tox_stc_observable_svd_workspace]].
        integer(int32), intent(in) :: iwork_size
            !! Size of tmp_iwork
            !! It is *VERY IMPORTANT* to compute this argument from the `iwork_size` output produced by [[tox_shape_truthful_clustering_observable_impl(module):tox_stc_observable_svd_workspace]].
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        logical(c_bool), dimension(n_vectors), intent(in) :: member_selection_mask
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

        call observable_impl(&
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
    end subroutine observable_expert

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_observable_impl(module):ensemble_final_observable_impl]].
    !| Not simply the last *populated* history column: `stc_push_ensemble_history`
    !| (`tox_shape_truthful_clustering_impl`) also pushes a *rejected* final candidate right
    !| before `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
    !| `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is, in exactly those
    !| two cases, the discarded candidate's geometry, not the ensemble's real final state.
    !| Scans each ensemble's history backward for the last column that is both populated
    !| (`ensemble_k_history /= 0`) and itself accepted (`ensemble_accepted_history`), and
    !| slices `U`/`d`/`S`/`mu`/`G`/`k` out at that column. `ensemble_has_final` is `.false.`
    !| (all other outputs zero for that ensemble) only when no column qualifies at all --
    !| possible for `STOP_REASON_MAX_SIZE` firing at the bootstrap step itself, before any
    !| genuine SVD ever ran for that seed, and rarely when a small `o` lets a rejected push
    !| evict every accepted entry the window ever held (`ensemble_U_first`/`ensemble_d_first`
    !| still hold the bootstrap iteration in that case, but this kernel does not fall back to
    !| them, since there is no `G_first`/`mu_first` counterpart to complete a fallback "final
    !| state" from -- see `misc/mod_STC.md`, "Ensemble identification"). Consolidates what was
    !| previously a private, single-consumer helper in `tox_stc_json.F90`
    !| (`stc_last_accepted_history_index`) into a proper, independently testable kernel now
    !| that `tox_shape_truthful_clustering_filter_impl` needs the exact same extraction --
    !| placed here, not in the parent module, specifically so both of those (and any future
    !| sibling) can depend on it without a circular module dependency (the parent module
    !| already `use`s `tox_shape_truthful_clustering_reconciliation_impl`, which will in turn
    !| `use` the filter kernel module, which needs this).
    pure subroutine ensemble_final_observable(&
            n_dimensions,&
            o,&
            n_ensembles,&
            ensemble_U_history,&
            ensemble_d_history,&
            ensemble_S_history,&
            ensemble_mu_history,&
            ensemble_G_history,&
            ensemble_k_history,&
            ensemble_accepted_history,&
            ensemble_U_final,&
            ensemble_d_final,&
            ensemble_S_final,&
            ensemble_mu_final,&
            ensemble_G_final,&
            ensemble_k_final,&
            ensemble_has_final,&
            ensemble_final_index,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_dimensions, n_dimensions, o, n_ensembles), intent(in) :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
            !! merged `ensemble_U_history`
        integer(int32), dimension(o, n_ensembles), intent(in) :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), dimension(n_dimensions, o, n_ensembles), intent(in) :: ensemble_S_history
            !! Per-ensemble trailing singular values
        real(real64), dimension(n_dimensions, o, n_ensembles), intent(in) :: ensemble_mu_history
            !! Per-ensemble trailing centers
        real(real64), dimension(o, n_ensembles), intent(in) :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        integer(int32), dimension(o, n_ensembles), intent(in) :: ensemble_k_history
            !! Per-ensemble trailing sizes; 0 marks an unpopulated column
        logical(c_bool), dimension(o, n_ensembles), intent(in) :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself
            !! accepted -- see this kernel's own summary above
        real(real64), dimension(n_dimensions, n_dimensions, n_ensembles), intent(out) :: ensemble_U_final
            !! Each ensemble's final accepted tangent+normal basis; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        integer(int32), dimension(n_ensembles), intent(out) :: ensemble_d_final
            !! Each ensemble's final accepted intrinsic dimension; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(real64), dimension(n_dimensions, n_ensembles), intent(out) :: ensemble_S_final
            !! Each ensemble's final accepted singular values; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(real64), dimension(n_dimensions, n_ensembles), intent(out) :: ensemble_mu_final
            !! Each ensemble's final accepted center; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        real(real64), dimension(n_ensembles), intent(out) :: ensemble_G_final
            !! Each ensemble's final accepted spectral gap; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        integer(int32), dimension(n_ensembles), intent(out) :: ensemble_k_final
            !! Each ensemble's final accepted size; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        logical(c_bool), dimension(n_ensembles), intent(out) :: ensemble_has_final
            !! Whether any history column at all qualifies as this ensemble's final accepted
            !! state -- see this kernel's own summary above for the (rare) `.false.` cases
        integer(int32), dimension(n_ensembles), intent(out) :: ensemble_final_index
            !! The history column each `_final` output was sliced from (0 when
            !! `ensemble_has_final` is `.false.`) -- also, since every column 1..this index is
            !! itself guaranteed accepted (only ever the single *last* populated column can be
            !! the rejected candidate this kernel's own summary describes), this doubles as the
            !! count of genuinely accepted, plottable history columns, for callers (e.g.
            !! `tox_stc_json`'s own `observable_history`) that need to iterate the whole
            !! trailing window, not just its final entry
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=1_int32, min=2_int32)
        call validate_in_range_int(o, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_ensembles, ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_real(ensemble_U_history, n_dimensions * n_dimensions * o * n_ensembles, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(ensemble_S_history, n_dimensions * o * n_ensembles, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(ensemble_mu_history, n_dimensions * o * n_ensembles, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(ensemble_G_history, o * n_ensembles, ierr, arg_pos=8_int32)
        if (is_err(ierr)) return
#endif

        call ensemble_final_observable_impl(&
            n_dimensions = n_dimensions,&
            o = o,&
            n_ensembles = n_ensembles,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history,&
            ensemble_U_final = ensemble_U_final,&
            ensemble_d_final = ensemble_d_final,&
            ensemble_S_final = ensemble_S_final,&
            ensemble_mu_final = ensemble_mu_final,&
            ensemble_G_final = ensemble_G_final,&
            ensemble_k_final = ensemble_k_final,&
            ensemble_has_final = ensemble_has_final,&
            ensemble_final_index = ensemble_final_index&
        )
    end subroutine ensemble_final_observable

end module tox_shape_truthful_clustering_observable
