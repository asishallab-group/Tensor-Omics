#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_observable(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_observable_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: normal_error_c
    public :: tangent_scales_c
    public :: observable_expert_c
    public :: observable_c
    public :: ensemble_final_observable_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable(module):normal_error(subroutine)]]
    !| No pass over the ensemble's member vectors is required; the sum is already implied by
    !| the singular value decomposition [[tox_shape_truthful_clustering_observable_kernel(module):observable_kernel]]
    !| computes.
    subroutine normal_error_c(&
            d,&
            eigenvalues,&
            n_dimensions,&
            normal_error_value,&
            ierr&
        ) bind(C, name="normal_error_c")
        use tox_shape_truthful_clustering_observable, only: normal_error

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: eigenvalues
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(out), target :: normal_error_value
            !! Mean squared residual off the d-dimensional tangent subspace
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(d)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(normal_error_value)
        M_CHECK_ARRAY_NON_NULL(eigenvalues, n_dimensions)

        call normal_error(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            normal_error_value = normal_error_value,&
            ierr = ierr&
        )
    end subroutine normal_error_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable(module):tangent_scales(subroutine)]]
    subroutine tangent_scales_c(&
            d,&
            eigenvalues,&
            n_dimensions,&
            tangent_scales_value,&
            ierr&
        ) bind(C, name="tangent_scales_c")
        use tox_shape_truthful_clustering_observable, only: tangent_scales

        integer(c_int), intent(in), target :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        real(c_double), dimension(n_dimensions), intent(in), target :: eigenvalues
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(d), intent(out), target :: tangent_scales_value
            !! Extent along each of the d tangent directions
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(d)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_ARRAY_NON_NULL(eigenvalues, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tangent_scales_value, d)

        call tangent_scales(&
            d = d,&
            eigenvalues = eigenvalues,&
            n_dimensions = n_dimensions,&
            tangent_scales_value = tangent_scales_value,&
            ierr = ierr&
        )
    end subroutine tangent_scales_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable(module):observable(subroutine)]]
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    subroutine observable_expert_c(&
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
        ) bind(C, name="observable_expert_c")
        use tox_shape_truthful_clustering_observable, only: observable

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_member
            !! Number of selected members (count of .TRUE. in member_selection_mask)
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: lwork
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork` output produced by [[tox_shape_truthful_clustering_observable_kernel(module):tox_stc_observable_svd_workspace]].
        integer(c_int), intent(in), target :: iwork_size
            !! Size of tmp_iwork
            !! It is *VERY IMPORTANT* to compute this argument from the `iwork_size` output produced by [[tox_shape_truthful_clustering_observable_kernel(module):tox_stc_observable_svd_workspace]].
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        logical(c_bool), dimension(n_vectors), intent(in), target :: member_selection_mask
            !! Ensemble membership over the full dataset
        real(c_double), dimension(n_dimensions, n_selected_member), intent(out), target :: tmp_y
            !! Workspace: centered member matrix
        real(c_double), dimension(min(n_dimensions,n_selected_member)), intent(out), target :: tmp_s
            !! Workspace: singular values
        real(c_double), dimension(n_dimensions, min(n_dimensions,n_selected_member)), intent(out), target :: tmp_u_econ
            !! Workspace: economy-mode left singular vectors
        real(c_double), dimension(min(n_dimensions,n_selected_member), n_selected_member), intent(out), target :: tmp_vt_econ
            !! Workspace: economy-mode right singular vectors, transposed (unused beyond the SVD call)
        real(c_double), dimension(lwork), intent(out), target :: tmp_work
            !! Workspace: LAPACK dgesdd scratch
        integer(c_int), dimension(iwork_size), intent(out), target :: tmp_iwork
            !! Workspace: LAPACK dgesdd integer scratch
        real(c_double), dimension(n_dimensions, n_dimensions), intent(out), target :: U
            !! Tangent+normal basis, zero-padded beyond rank
        real(c_double), dimension(n_dimensions), intent(out), target :: eigenvalues
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(c_double), dimension(n_dimensions), intent(out), target :: mu
            !! Ensemble center
        integer(c_int), intent(out), target :: d
            !! Estimated intrinsic (tangent) dimension
        real(c_double), intent(out), target :: G
            !! Spectral gap at d
        real(c_double), intent(out), target :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(c_double), dimension(n_dimensions), intent(out), target :: tangent_scales_value
            !! Extent along each tangent direction, zero-padded beyond d
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        logical, dimension(n_vectors) :: member_selection_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_member)
        M_CHECK_NON_NULL(lwork)
        M_CHECK_NON_NULL(iwork_size)
        M_CHECK_NON_NULL(d)
        M_CHECK_NON_NULL(G)
        M_CHECK_NON_NULL(normal_error_value)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(member_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_y, n_dimensions * n_selected_member)
        M_CHECK_ARRAY_NON_NULL(tmp_s, (min(n_dimensions,n_selected_member)))
        M_CHECK_ARRAY_NON_NULL(tmp_u_econ, n_dimensions * (min(n_dimensions,n_selected_member)))
        M_CHECK_ARRAY_NON_NULL(tmp_vt_econ, (min(n_dimensions,n_selected_member)) * n_selected_member)
        M_CHECK_ARRAY_NON_NULL(tmp_work, lwork)
        M_CHECK_ARRAY_NON_NULL(tmp_iwork, iwork_size)
        M_CHECK_ARRAY_NON_NULL(U, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(eigenvalues, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(mu, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tangent_scales_value, n_dimensions)

        member_selection_mask_f = member_selection_mask

        call observable(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            member_selection_mask = member_selection_mask_f,&
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
    end subroutine observable_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable(module):observable_alloc(subroutine)]]
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    subroutine observable_c(&
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
        ) bind(C, name="observable_c")
        use tox_shape_truthful_clustering_observable, only: observable_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        logical(c_bool), dimension(n_vectors), intent(in), target :: member_selection_mask
            !! Ensemble membership over the full dataset
        integer(c_int), intent(in), target :: n_selected_member
            !! Number of selected members (count of .TRUE. in member_selection_mask)
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_dimensions, n_dimensions), intent(out), target :: U
            !! Tangent+normal basis, zero-padded beyond rank
        real(c_double), dimension(n_dimensions), intent(out), target :: eigenvalues
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(c_double), dimension(n_dimensions), intent(out), target :: mu
            !! Ensemble center
        integer(c_int), intent(out), target :: d
            !! Estimated intrinsic (tangent) dimension
        real(c_double), intent(out), target :: G
            !! Spectral gap at d
        real(c_double), intent(out), target :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(c_double), dimension(n_dimensions), intent(out), target :: tangent_scales_value
            !! Extent along each tangent direction, zero-padded beyond d
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        logical, dimension(n_vectors) :: member_selection_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_member)
        M_CHECK_NON_NULL(d)
        M_CHECK_NON_NULL(G)
        M_CHECK_NON_NULL(normal_error_value)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(member_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(U, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(eigenvalues, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(mu, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tangent_scales_value, n_dimensions)

        member_selection_mask_f = member_selection_mask

        call observable_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            member_selection_mask = member_selection_mask_f,&
            n_selected_member = n_selected_member,&
            U = U,&
            eigenvalues = eigenvalues,&
            mu = mu,&
            d = d,&
            G = G,&
            normal_error_value = normal_error_value,&
            tangent_scales_value = tangent_scales_value,&
            ierr = ierr&
        )
    end subroutine observable_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_observable(module):ensemble_final_observable(subroutine)]]
    !| Not simply the last *populated* history column: `stc_push_ensemble_history`
    !| (`tox_shape_truthful_clustering_kernel`) also pushes a *rejected* final candidate right
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
    !| that `tox_shape_truthful_clustering_filter_kernel` needs the exact same extraction --
    !| placed here, not in the parent module, specifically so both of those (and any future
    !| sibling) can depend on it without a circular module dependency (the parent module
    !| already `use`s `tox_shape_truthful_clustering_reconciliation_kernel`, which will in turn
    !| `use` the filter kernel module, which needs this).
    subroutine ensemble_final_observable_c(&
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
        ) bind(C, name="ensemble_final_observable_c")
        use tox_shape_truthful_clustering_observable, only: ensemble_final_observable

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_dimensions, n_dimensions, o, n_ensembles), intent(in), target :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
            !! merged `ensemble_U_history`
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_S_history
            !! Per-ensemble trailing singular values
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_mu_history
            !! Per-ensemble trailing centers
        real(c_double), dimension(o, n_ensembles), intent(in), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_k_history
            !! Per-ensemble trailing sizes; 0 marks an unpopulated column
        logical(c_bool), dimension(o, n_ensembles), intent(in), target :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself
            !! accepted -- see this kernel's own summary above
        real(c_double), dimension(n_dimensions, n_dimensions, n_ensembles), intent(out), target :: ensemble_U_final
            !! Each ensemble's final accepted tangent+normal basis; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        integer(c_int), dimension(n_ensembles), intent(out), target :: ensemble_d_final
            !! Each ensemble's final accepted intrinsic dimension; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(c_double), dimension(n_dimensions, n_ensembles), intent(out), target :: ensemble_S_final
            !! Each ensemble's final accepted singular values; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(c_double), dimension(n_dimensions, n_ensembles), intent(out), target :: ensemble_mu_final
            !! Each ensemble's final accepted center; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        real(c_double), dimension(n_ensembles), intent(out), target :: ensemble_G_final
            !! Each ensemble's final accepted spectral gap; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        integer(c_int), dimension(n_ensembles), intent(out), target :: ensemble_k_final
            !! Each ensemble's final accepted size; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        logical(c_bool), dimension(n_ensembles), intent(out), target :: ensemble_has_final
            !! Whether any history column at all qualifies as this ensemble's final accepted
            !! state -- see this kernel's own summary above for the (rare) `.false.` cases
        integer(c_int), dimension(n_ensembles), intent(out), target :: ensemble_final_index
            !! The history column each `_final` output was sliced from (0 when
            !! `ensemble_has_final` is `.false.`) -- also, since every column 1..this index is
            !! itself guaranteed accepted (only ever the single *last* populated column can be
            !! the rejected candidate this kernel's own summary describes), this doubles as the
            !! count of genuinely accepted, plottable history columns, for callers (e.g.
            !! `tox_stc_json`'s own `observable_history`) that need to iterate the whole
            !! trailing window, not just its final entry
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(o, n_ensembles) :: ensemble_accepted_history_f
        logical, dimension(n_ensembles) :: ensemble_has_final_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_history, n_dimensions * n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_history, n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_history, n_dimensions * o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_accepted_history, o * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_final, n_dimensions * n_dimensions * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_final, n_dimensions * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_final, n_dimensions * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_has_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_final_index, n_ensembles)

        ensemble_accepted_history_f = ensemble_accepted_history

        call ensemble_final_observable(&
            n_dimensions = n_dimensions,&
            o = o,&
            n_ensembles = n_ensembles,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history_f,&
            ensemble_U_final = ensemble_U_final,&
            ensemble_d_final = ensemble_d_final,&
            ensemble_S_final = ensemble_S_final,&
            ensemble_mu_final = ensemble_mu_final,&
            ensemble_G_final = ensemble_G_final,&
            ensemble_k_final = ensemble_k_final,&
            ensemble_has_final = ensemble_has_final_f,&
            ensemble_final_index = ensemble_final_index,&
            ierr = ierr&
        )

        ensemble_has_final = ensemble_has_final_f
    end subroutine ensemble_final_observable_c

end module tox_shape_truthful_clustering_observable_c
#endif
