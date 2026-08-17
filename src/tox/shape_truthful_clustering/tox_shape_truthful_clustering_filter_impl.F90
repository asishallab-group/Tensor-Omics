#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Ensemble Filtering
!|
!| `filter_ensembles`: decides which ensembles are eligible to be submitted to
!| `merge_to_super_ensembles` (see the sibling `tox_shape_truthful_clustering_reconciliation_impl`,
!| whose own `ensemble_reconciliation` is now a thin two-call orchestrator: this module first,
!| then that one). Composed of independent, individually testable per-criterion filters --
!| `filter_ensembles_by_stop_condition`, `filter_ensembles_by_dimension`,
!| `filter_ensembles_by_var_explained` -- each returning its own `eligible(n_ensembles)` mask
!| over the *same* ensembles, combined by `filter_ensembles` itself via a plain logical AND. A
!| criterion whose own threshold/allowed-set argument is omitted contributes an all-`.true.`
!| mask (no constraint from that criterion), so omitting every optional argument makes
!| `filter_ensembles` itself a true no-op (every ensemble eligible) -- see `misc/mod_STC.md`,
!| "Ensemble Reconciliation".
!|
!| Filtering never alters `ensemble_identification`'s own output, nor does it remove an
!| ineligible ensemble from anywhere else this whole family reports it (points, the JSON's
!| `ensembles` array, CSV output, ...) -- only `merge_to_super_ensembles`'s own pairing/grouping
!| decision (and, downstream, `tox_stc_json`'s independently-computed `overlap_coefficient_matrix`,
!| which applies the identical mask for the same reason) ever sees an ineligible ensemble
!| excluded. No array copying/compaction anywhere in this module: every mask is exactly
!| `n_ensembles` long, over the same 1-indexed ensemble numbering everything else in this
!| family already uses.
module tox_shape_truthful_clustering_filter_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_shape_truthful_clustering_observable_impl, only: ensemble_final_observable_impl, normal_error_impl
    M_IMPLICIT_NONE

    private
    public :: filter_ensembles_by_stop_condition_impl
    public :: filter_ensembles_by_dimension_impl
    public :: filter_ensembles_by_var_explained_impl
    public :: filter_ensembles_impl

contains

    !> summary: Ensemble eligibility by Stop Condition
    !| AUTHOR_ASIS_HALLAB
    !| `eligible(e) = allowed_stop_reasons(ensemble_stop_reason(e))` when `allowed_stop_reasons`
    !| is present; all `.true.` (no filtering) when absent.
    pure subroutine filter_ensembles_by_stop_condition_impl(n_ensembles, ensemble_stop_reason, &
                                                               allowed_stop_reasons, eligible)
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: ensemble_stop_reason(n_ensembles)
            !! Per-ensemble Stop Condition, see `ensemble_identification`'s merged
            !! `ensemble_stop_reason` -- an index 1..4 into `allowed_stop_reasons` below, in the
            !! order `tox_shape_truthful_clustering_impl`'s own `STOP_REASON_MAX_SIZE` (1),
            !! `STOP_REASON_REJECTED_AFTER_STABLE` (2), `STOP_REASON_REJECTED_IMMEDIATELY` (3),
            !! `STOP_REASON_FIXED_POINT` (4) -- not imported by name here, to avoid a circular
            !! module dependency (the parent module already `use`s the reconciliation module,
            !! which `use`s this one)
            !! DM_MIN(1_int32)
            !! DM_MAX(4_int32)
        logical(c_bool), intent(in), optional :: allowed_stop_reasons(4)
            !! Per-Stop-Condition eligibility, indexed as documented on `ensemble_stop_reason`
            !! above. Absent means no filtering (every Stop Condition allowed) -- deliberately
            !! nullable, not annotated with a generated default: the generator only evaluates
            !! constant *scalar* expressions for that annotation (`codegen_guide.md` section 5.5)
        logical(c_bool), intent(out) :: eligible(n_ensembles)
            !! Per-ensemble eligibility from this criterion alone

        integer(int32) :: e
        logical(c_bool)        :: actual_allowed_stop_reasons(4)

        if (present(allowed_stop_reasons)) then
            actual_allowed_stop_reasons = allowed_stop_reasons
        else
            actual_allowed_stop_reasons = .true.
        end if

        do e = 1, n_ensembles
            eligible(e) = actual_allowed_stop_reasons(ensemble_stop_reason(e))
        end do

    end subroutine filter_ensembles_by_stop_condition_impl

    !> summary: Ensemble eligibility by final intrinsic dimension
    !| AUTHOR_ASIS_HALLAB
    !| `d_min <= d <= d_max`, both inclusive, each independently optional (an absent bound
    !| contributes no constraint on that side). An ensemble with no final accepted state at all
    !| (`ensemble_has_final(e)` false) is never eligible under this criterion once at least one
    !| of `d_min`/`d_max` is supplied -- there is no `d` to judge. Both bounds absent is a true
    !| no-op (every ensemble eligible, `ensemble_has_final` not even consulted), matching
    !| `filter_ensembles_by_stop_condition_impl`'s own "omitted means unconstrained"
    !| convention. `d_min` could in principle default to `0_int32` (a genuine constant
    !| expression), but is left nullable like `d_max` (whose own natural default, `n_dimensions`,
    !| is a runtime value and so cannot be a generated default) rather than have the two bounds
    !| of one range behave asymmetrically.
    pure subroutine filter_ensembles_by_dimension_impl(n_dimensions, n_ensembles, ensemble_d_final, &
                                                          ensemble_has_final, d_min, d_max, eligible)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: ensemble_d_final(n_ensembles)
            !! Each ensemble's final accepted intrinsic dimension, see
            !! `ensemble_final_observable`
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        logical(c_bool), intent(in) :: ensemble_has_final(n_ensembles)
            !! Whether each ensemble has a final accepted state at all, see
            !! `ensemble_final_observable`
        integer(int32), intent(in), optional :: d_min
            !! Minimum tolerated final intrinsic dimension, inclusive
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in), optional :: d_max
            !! Maximum tolerated final intrinsic dimension, inclusive
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        logical(c_bool), intent(out) :: eligible(n_ensembles)
            !! Per-ensemble eligibility from this criterion alone

        integer(int32) :: e

        if (.not. present(d_min) .and. .not. present(d_max)) then
            eligible = .true.
            return
        end if

        do e = 1, n_ensembles
            eligible(e) = ensemble_has_final(e)
            if (.not. eligible(e)) cycle
            if (present(d_min)) eligible(e) = eligible(e) .and. (ensemble_d_final(e) >= d_min)
            if (present(d_max)) eligible(e) = eligible(e) .and. (ensemble_d_final(e) <= d_max)
        end do

    end subroutine filter_ensembles_by_dimension_impl

    !> summary: Ensemble eligibility by final classical variance explained
    !| AUTHOR_ASIS_HALLAB
    !| Variance explained $= \sum_{j=1}^{d}\lambda_j / (\sum_{j=1}^{d}\lambda_j +
    !| \text{normal\_error})$, the familiar PCA energy-ratio (not the scale/sqrt-based ratio
    !| considered and rejected during this module's own design -- the classical, squared-units
    !| form was chosen specifically for being the measure data scientists already expect), with
    !| eigenvalues recovered from the final singular values exactly as `observable` itself does,
    !| $\lambda_j = s_j^2/(k-1)$, and `normal_error` reusing
    !| [[tox_shape_truthful_clustering_observable_impl(module):normal_error_impl]] directly
    !| rather than re-deriving its sum. An ensemble with no final accepted state, or whose final
    !| size is too small for the $k-1$ denominator to be meaningful ($k \leq 1$), is never
    !| eligible under this criterion once `var_explained_min` is supplied. Absent
    !| `var_explained_min` is a true no-op, matching this module's other two filters.
    pure subroutine filter_ensembles_by_var_explained_impl(n_dimensions, n_ensembles, ensemble_S_final, &
                                                              ensemble_d_final, ensemble_k_final, &
                                                              ensemble_has_final, var_explained_min, eligible)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_ensembles
            !! Number of ensembles N_E
            !! DM_MIN(0_int32)
        real(real64), intent(in) :: ensemble_S_final(n_dimensions, n_ensembles)
            !! Each ensemble's final accepted singular values, see `ensemble_final_observable`
        integer(int32), intent(in) :: ensemble_d_final(n_ensembles)
            !! Each ensemble's final accepted intrinsic dimension
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: ensemble_k_final(n_ensembles)
            !! Each ensemble's final accepted size
            !! DM_MIN(0_int32)
        logical(c_bool), intent(in) :: ensemble_has_final(n_ensembles)
            !! Whether each ensemble has a final accepted state at all
        real(real64), intent(in), optional :: var_explained_min
            !! Minimum tolerated fraction of variance explained by the tangent subspace,
            !! inclusive
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: eligible(n_ensembles)
            !! Per-ensemble eligibility from this criterion alone

        integer(int32) :: e, d_e
        real(real64)   :: eigenvalues(n_dimensions), tangent_sum, normal_error_value, variance_explained

        if (.not. present(var_explained_min)) then
            eligible = .true.
            return
        end if

        do e = 1, n_ensembles
            if (.not. ensemble_has_final(e) .or. ensemble_k_final(e) <= 1) then
                eligible(e) = .false.
                cycle
            end if

            d_e = ensemble_d_final(e)
            eigenvalues = ensemble_S_final(:, e)**2/real(ensemble_k_final(e) - 1, real64)
            tangent_sum = sum(eigenvalues(1:d_e))
            call normal_error_impl(d_e, eigenvalues, n_dimensions, normal_error_value)
            variance_explained = tangent_sum/(tangent_sum + normal_error_value + epsilon(1.0_real64))
            eligible(e) = variance_explained >= var_explained_min
        end do

    end subroutine filter_ensembles_by_var_explained_impl

    !> summary: Combined ensemble eligibility for `merge_to_super_ensembles`
    !| AUTHOR_ASIS_HALLAB
    !| Extracts each ensemble's final accepted state once (via `ensemble_final_observable`,
    !| shared so this and `tox_stc_json`'s own reporting never derive it two different ways),
    !| then calls each of this module's three per-criterion filters and combines their masks
    !| with a plain logical AND. Also returns the three individual masks, not just the
    !| combination -- so a caller (the report, in particular) can say *which* criterion excluded
    !| a given ensemble, not merely that one did. Supplying none of `allowed_stop_reasons`/
    !| `d_min`/`d_max`/`var_explained_min` makes every ensemble eligible (all four masks
    !| all-`.true.`), matching each individual filter's own no-op convention.
    pure subroutine filter_ensembles_impl(n_dimensions, o, n_ensembles, &
                                            ensemble_U_history, ensemble_d_history, ensemble_S_history, &
                                            ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                                            ensemble_accepted_history, ensemble_stop_reason, &
                                            allowed_stop_reasons, d_min, d_max, var_explained_min, &
                                            eligible, eligible_by_stop_condition, eligible_by_dimension, &
                                            eligible_by_var_explained)
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
            !! merged output
        integer(int32), intent(in) :: ensemble_d_history(o, n_ensembles)
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), intent(in) :: ensemble_S_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing singular values
        real(real64), intent(in) :: ensemble_mu_history(n_dimensions, o, n_ensembles)
            !! Per-ensemble trailing centers
        real(real64), intent(in) :: ensemble_G_history(o, n_ensembles)
            !! Per-ensemble trailing spectral gaps
        integer(int32), intent(in) :: ensemble_k_history(o, n_ensembles)
            !! Per-ensemble trailing sizes
        logical(c_bool), intent(in) :: ensemble_accepted_history(o, n_ensembles)
            !! Whether the growth iteration retained in each history column was itself accepted
        integer(int32), intent(in) :: ensemble_stop_reason(n_ensembles)
            !! Per-ensemble Stop Condition, see `filter_ensembles_by_stop_condition_impl`
            !! DM_MIN(1_int32)
            !! DM_MAX(4_int32)
        logical(c_bool), intent(in), optional :: allowed_stop_reasons(4)
            !! See `filter_ensembles_by_stop_condition_impl`
        integer(int32), intent(in), optional :: d_min
            !! See `filter_ensembles_by_dimension_impl`
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in), optional :: d_max
            !! See `filter_ensembles_by_dimension_impl`
            !! DM_MIN(0_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in), optional :: var_explained_min
            !! See `filter_ensembles_by_var_explained_impl`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: eligible(n_ensembles)
            !! Combined eligibility: `.true.` only where all three per-criterion masks are
        logical(c_bool), intent(out) :: eligible_by_stop_condition(n_ensembles)
            !! See `filter_ensembles_by_stop_condition_impl`
        logical(c_bool), intent(out) :: eligible_by_dimension(n_ensembles)
            !! See `filter_ensembles_by_dimension_impl`
        logical(c_bool), intent(out) :: eligible_by_var_explained(n_ensembles)
            !! See `filter_ensembles_by_var_explained_impl`

        real(real64)   :: ensemble_U_final(n_dimensions, n_dimensions, n_ensembles)
        integer(int32) :: ensemble_d_final(n_ensembles)
        real(real64)   :: ensemble_S_final(n_dimensions, n_ensembles)
        real(real64)   :: ensemble_mu_final(n_dimensions, n_ensembles)
        real(real64)   :: ensemble_G_final(n_ensembles)
        integer(int32) :: ensemble_k_final(n_ensembles)
        logical(c_bool)        :: ensemble_has_final(n_ensembles)
        integer(int32) :: ensemble_final_index(n_ensembles)

        call ensemble_final_observable_impl(n_dimensions, o, n_ensembles, &
            ensemble_U_history, ensemble_d_history, ensemble_S_history, &
            ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, &
            ensemble_U_final, ensemble_d_final, ensemble_S_final, ensemble_mu_final, ensemble_G_final, &
            ensemble_k_final, ensemble_has_final, ensemble_final_index)

        call filter_ensembles_by_stop_condition_impl(n_ensembles, ensemble_stop_reason, allowed_stop_reasons, &
            eligible_by_stop_condition)

        call filter_ensembles_by_dimension_impl(n_dimensions, n_ensembles, ensemble_d_final, ensemble_has_final, &
            d_min, d_max, eligible_by_dimension)

        call filter_ensembles_by_var_explained_impl(n_dimensions, n_ensembles, ensemble_S_final, ensemble_d_final, &
            ensemble_k_final, ensemble_has_final, var_explained_min, eligible_by_var_explained)

        eligible = eligible_by_stop_condition .and. eligible_by_dimension .and. eligible_by_var_explained

    end subroutine filter_ensembles_impl

end module tox_shape_truthful_clustering_filter_impl
