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
module tox_shape_truthful_clustering_observable_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, ERR_INTERNAL
    M_IMPLICIT_NONE

    interface
        ! Declared pure to document that dgesdd is a deterministic numerical routine with no
        ! I/O and no global state -- which is what makes it safe to call from a pure kernel;
        ! not a language requirement, just documentation of its own thread-safety (matching
        ! dsyev's identical precedent in src/lomanle.F90).
        pure subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            import :: int32, real64
            character,      intent(in)    :: jobz
            integer(int32), intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(real64),   intent(inout) :: a(lda, n)
            real(real64),   intent(out)   :: s(min(m, n))
            real(real64),   intent(out)   :: u(ldu, min(m, n))
            real(real64),   intent(out)   :: vt(ldvt, n)
            real(real64),   intent(out)   :: work(lwork)
            integer(int32), intent(out)   :: iwork(8*min(m, n))
            integer(int32), intent(out)   :: info
        end subroutine dgesdd
    end interface

    private
    public :: normal_error_impl
    public :: tangent_scales_impl
    public :: observable_impl
    public :: tox_stc_observable_svd_workspace
    public :: ensemble_final_observable_impl

contains

    !> summary: Mean squared residual of an ensemble's members off its tangent subspace
    !| AUTHOR_ASIS_HALLAB
    !| No pass over the ensemble's member vectors is required; the sum is already implied by
    !| the singular value decomposition [[tox_shape_truthful_clustering_observable_impl(module):observable_impl]]
    !| computes.
    pure subroutine normal_error_impl(d, eigenvalues, n_dimensions, normal_error_value)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! DM_MIN(0.0_real64)
        real(real64), intent(out) :: normal_error_value
            !! Mean squared residual off the d-dimensional tangent subspace

        normal_error_value = sum(eigenvalues(d + 1:n_dimensions))

    end subroutine normal_error_impl

    !> summary: Extent along each tangent direction of an ensemble's tangent subspace
    !| AUTHOR_ASIS_HALLAB
    pure subroutine tangent_scales_impl(d, eigenvalues, n_dimensions, tangent_scales_value)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: d
            !! Intrinsic (tangent) dimension of the ensemble
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: eigenvalues(n_dimensions)
            !! Ensemble covariance eigenvalues, descending: lambda_1 >= ... >= lambda_D >= 0
            !! DM_MIN(0.0_real64)
        real(real64), intent(out) :: tangent_scales_value(d)
            !! Extent along each of the d tangent directions

        tangent_scales_value = sqrt(eigenvalues(1:d))

    end subroutine tangent_scales_impl

    !> M_EXPORT_C
    !| summary: Recommend LAPACK dgesdd workspace sizes for observable's economy-mode SVD
    !| AUTHOR_ASIS_HALLAB
    !| The documented minimum-workspace formula for JOBZ='S' (see `man dgesdd`):
    !| LWORK >= 4*min(M,N)**2 + 7*min(M,N), IWORK size = 8*min(M,N), where M=n_dimensions and
    !| N=n_selected_member.
    pure subroutine tox_stc_observable_svd_workspace(n_dimensions, n_selected_member, lwork, iwork_size)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_selected_member
            !! Number of selected ensemble members
        integer(int32), intent(out) :: lwork
            !! Recommended size of the real LAPACK workspace
        integer(int32), intent(out) :: iwork_size
            !! Recommended size of the integer LAPACK workspace

        integer(int32) :: rank

        rank = min(n_dimensions, n_selected_member)
        lwork = 4*rank*rank + 7*rank
        iwork_size = 8*rank

    end subroutine tox_stc_observable_svd_workspace

    !> summary: The tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble
    !| AUTHOR_ASIS_HALLAB
    !| `U` and `eigenvalues` are zero-padded to the full ambient dimension `n_dimensions`:
    !| the economy SVD only yields `rank = min(n_dimensions, n_selected_member)` genuine
    !| columns/values, less than `n_dimensions` whenever an ensemble is smaller than the
    !| ambient space (typical early in growth). This keeps the output shape fixed regardless
    !| of ensemble size, and slots directly into `normal_error`/`tangent_scales`'s existing
    !| `n_dimensions`-length interface. `ierr` is set only if the LAPACK SVD fails to
    !| converge -- not a condition any input check could foresee.
    pure subroutine observable_impl(vectors, n_dimensions, n_vectors, member_selection_mask, n_selected_member, &
                                      lwork, iwork_size, &
                                      tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                                      U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        logical(c_bool), intent(in) :: member_selection_mask(n_vectors)
            !! Ensemble membership over the full dataset
        integer(int32), intent(in) :: n_selected_member
            !! Number of selected members (count of .TRUE. in member_selection_mask)
            !! DM_MIN(2_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: lwork
            !! Size of tmp_work
            !! DM_OUTPUT_FROM(lwork, tox_stc_observable_svd_workspace, tox_shape_truthful_clustering_observable_impl, AUTO)
        integer(int32), intent(in) :: iwork_size
            !! Size of tmp_iwork
            !! DM_OUTPUT_FROM(iwork_size, tox_stc_observable_svd_workspace, tox_shape_truthful_clustering_observable_impl, AUTO)
        real(real64), intent(out) :: tmp_y(n_dimensions, n_selected_member)
            !! Workspace: centered member matrix
        real(real64), intent(out) :: tmp_s(min(n_dimensions, n_selected_member))
            !! Workspace: singular values
        real(real64), intent(out) :: tmp_u_econ(n_dimensions, min(n_dimensions, n_selected_member))
            !! Workspace: economy-mode left singular vectors
        real(real64), intent(out) :: tmp_vt_econ(min(n_dimensions, n_selected_member), n_selected_member)
            !! Workspace: economy-mode right singular vectors, transposed (unused beyond the SVD call)
        real(real64), intent(out) :: tmp_work(lwork)
            !! Workspace: LAPACK dgesdd scratch
        integer(int32), intent(out) :: tmp_iwork(iwork_size)
            !! Workspace: LAPACK dgesdd integer scratch
        real(real64), intent(out) :: U(n_dimensions, n_dimensions)
            !! Tangent+normal basis, zero-padded beyond rank
        real(real64), intent(out) :: eigenvalues(n_dimensions)
            !! Covariance eigenvalues, descending, zero-padded beyond rank
        real(real64), intent(out) :: mu(n_dimensions)
            !! Ensemble center
        integer(int32), intent(out) :: d
            !! Estimated intrinsic (tangent) dimension
        real(real64), intent(out) :: G
            !! Spectral gap at d
        real(real64), intent(out) :: normal_error_value
            !! Mean squared residual off the tangent subspace
        real(real64), intent(out) :: tangent_scales_value(n_dimensions)
            !! Extent along each tangent direction, zero-padded beyond d
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: rank, i_vec, i_member, j, r, info, d_local
        real(real64)   :: gap, best_gap

        call set_ok(ierr)

        rank = min(n_dimensions, n_selected_member)

        i_member = 0
        do i_vec = 1, n_vectors
            if (.not. member_selection_mask(i_vec)) cycle
            i_member = i_member + 1
            tmp_y(:, i_member) = vectors(:, i_vec)
        end do

        mu    = sum(tmp_y, dim=2)/real(n_selected_member, real64)
        tmp_y = tmp_y - spread(mu, dim=2, ncopies=n_selected_member)

        call dgesdd('S', n_dimensions, n_selected_member, tmp_y, n_dimensions, tmp_s, &
                   tmp_u_econ, n_dimensions, tmp_vt_econ, rank, tmp_work, lwork, tmp_iwork, info)
        if (info /= 0) then
            call set_err_once(ierr, ERR_INTERNAL)
            return
        end if

        U = 0.0_real64
        U(:, 1:rank) = tmp_u_econ(:, 1:rank)

        eigenvalues = 0.0_real64
        do j = 1, rank
            eigenvalues(j) = tmp_s(j)**2/real(n_selected_member - 1, real64)
        end do

        ! Spectral gap G(r) = lambda_r / (lambda_{r+1} + eps), r = 1..D-1; d = argmax_r G(r).
        ! eigenvalues is zero-padded beyond `rank`, so any r >= rank+1 contributes
        ! G(r) = 0/eps = 0 -- never a spurious maximum -- while r = rank itself, whenever
        ! rank < n_dimensions, divides a genuine nonzero eigenvalue by (0 + eps) and
        ! correctly wins whenever the ensemble's data does not span the full ambient space
        ! (e.g. still few members, early in growth).
        best_gap = -1.0_real64
        d_local  = 1
        do r = 1, n_dimensions - 1
            gap = eigenvalues(r)/(eigenvalues(r + 1) + epsilon(1.0_real64))
            if (gap > best_gap) then
                best_gap = gap
                d_local  = r
            end if
        end do
        d = d_local
        G = best_gap

        call normal_error_impl(d_local, eigenvalues, n_dimensions, normal_error_value)

        tangent_scales_value = 0.0_real64
        block
            real(real64) :: tangent_scales_local(d_local)
            call tangent_scales_impl(d_local, eigenvalues, n_dimensions, tangent_scales_local)
            tangent_scales_value(1:d_local) = tangent_scales_local
        end block

    end subroutine observable_impl

    !> summary: Each ensemble's final *accepted* growth-history state
    !| AUTHOR_ASIS_HALLAB
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
    pure subroutine ensemble_final_observable_impl(n_dimensions, o, n_ensembles, &
                                                      ensemble_U_history, ensemble_d_history, ensemble_S_history, &
                                                      ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                                                      ensemble_accepted_history, &
                                                      ensemble_U_final, ensemble_d_final, ensemble_S_final, &
                                                      ensemble_mu_final, ensemble_G_final, ensemble_k_final, &
                                                      ensemble_has_final, ensemble_final_index)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: o
            !! Trailing observable-history window depth
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: ensemble_U_history(n_dimensions, n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing tangent+normal bases, see `ensemble_identification`'s
            !! merged `ensemble_U_history`
        integer(int32), intent(in) :: ensemble_d_history(o, n_ensembles)
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), intent(in) :: ensemble_S_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing singular values
        real(real64), intent(in) :: ensemble_mu_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing centers
        real(real64), intent(in) :: ensemble_G_history(o, n_ensembles)
            !! Per-ensemble trailing spectral gaps
        integer(int32), intent(in) :: ensemble_k_history(o, n_ensembles)
            !! Per-ensemble trailing sizes; 0 marks an unpopulated column
        logical(c_bool), intent(in) :: ensemble_accepted_history(o, n_ensembles)
            !! Whether the growth iteration retained in each history column was itself
            !! accepted -- see this kernel's own summary above
        real(real64), intent(out) :: ensemble_U_final(n_dimensions, n_dimensions, n_ensembles)
            !! Each ensemble's final accepted tangent+normal basis; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        integer(int32), intent(out) :: ensemble_d_final(n_ensembles)
            !! Each ensemble's final accepted intrinsic dimension; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(real64), intent(out) :: ensemble_S_final(n_dimensions, n_ensembles)
            !! Each ensemble's final accepted singular values; zero when
            !! `ensemble_has_final` is `.false.` for that ensemble
        real(real64), intent(out) :: ensemble_mu_final(n_dimensions, n_ensembles)
            !! Each ensemble's final accepted center; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        real(real64), intent(out) :: ensemble_G_final(n_ensembles)
            !! Each ensemble's final accepted spectral gap; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        integer(int32), intent(out) :: ensemble_k_final(n_ensembles)
            !! Each ensemble's final accepted size; zero when `ensemble_has_final` is
            !! `.false.` for that ensemble
        logical(c_bool), intent(out) :: ensemble_has_final(n_ensembles)
            !! Whether any history column at all qualifies as this ensemble's final accepted
            !! state -- see this kernel's own summary above for the (rare) `.false.` cases
        integer(int32), intent(out) :: ensemble_final_index(n_ensembles)
            !! The history column each `_final` output was sliced from (0 when
            !! `ensemble_has_final` is `.false.`) -- also, since every column 1..this index is
            !! itself guaranteed accepted (only ever the single *last* populated column can be
            !! the rejected candidate this kernel's own summary describes), this doubles as the
            !! count of genuinely accepted, plottable history columns, for callers (e.g.
            !! `tox_stc_json`'s own `observable_history`) that need to iterate the whole
            !! trailing window, not just its final entry

        integer(int32) :: e, j, idx

        do e = 1, n_ensembles
            idx = 0
            do j = o, 1, -1
                if (ensemble_k_history(j, e) /= 0 .and. ensemble_accepted_history(j, e)) then
                    idx = j
                    exit
                end if
            end do

            ensemble_final_index(e) = idx
            ensemble_has_final(e) = idx > 0
            if (idx > 0) then
                ensemble_U_final(:, :, e)  = ensemble_U_history(:, :, idx, e)
                ensemble_d_final(e)        = ensemble_d_history(idx, e)
                ensemble_S_final(:, e)     = ensemble_S_history(:, idx, e)
                ensemble_mu_final(:, e)    = ensemble_mu_history(:, idx, e)
                ensemble_G_final(e)        = ensemble_G_history(idx, e)
                ensemble_k_final(e)        = ensemble_k_history(idx, e)
            else
                ensemble_U_final(:, :, e)  = 0.0_real64
                ensemble_d_final(e)        = 0
                ensemble_S_final(:, e)     = 0.0_real64
                ensemble_mu_final(:, e)    = 0.0_real64
                ensemble_G_final(e)        = 0.0_real64
                ensemble_k_final(e)        = 0
            end if
        end do

    end subroutine ensemble_final_observable_impl

end module tox_shape_truthful_clustering_observable_impl
