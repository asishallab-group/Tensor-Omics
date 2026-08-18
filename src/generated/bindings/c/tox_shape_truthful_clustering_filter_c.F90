#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_filter(module)]]
!| # Shape Truthful Clustering (STC): Ensemble Filtering
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
module tox_shape_truthful_clustering_filter_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: filter_ensembles_by_stop_condition_c
    public :: filter_ensembles_by_dimension_c
    public :: filter_ensembles_by_var_explained_c
    public :: filter_ensembles_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_filter(module):filter_ensembles_by_stop_condition(subroutine)]]
    !| `eligible(e) = allowed_stop_reasons(ensemble_stop_reason(e))` when `allowed_stop_reasons`
    !| is present; all `.true.` (no filtering) when absent.
    subroutine filter_ensembles_by_stop_condition_c(&
            n_ensembles,&
            ensemble_stop_reason,&
            allowed_stop_reasons,&
            eligible,&
            ierr&
        ) bind(C, name="filter_ensembles_by_stop_condition_c")
        use tox_shape_truthful_clustering_filter, only: filter_ensembles_by_stop_condition

        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `ensemble_identification`'s merged
            !! `ensemble_stop_reason` -- an index 1..4 into `allowed_stop_reasons` below, in the
            !! order `tox_shape_truthful_clustering_impl`'s own `STOP_REASON_MAX_SIZE` (1),
            !! `STOP_REASON_REJECTED_AFTER_STABLE` (2), `STOP_REASON_REJECTED_IMMEDIATELY` (3),
            !! `STOP_REASON_FIXED_POINT` (4) -- not imported by name here, to avoid a circular
            !! module dependency (the parent module already `use`s the reconciliation module,
            !! which `use`s this one)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `4_int32`.
        logical(c_bool), dimension(4), intent(in), optional :: allowed_stop_reasons
            !! Per-Stop-Condition eligibility, indexed as documented on `ensemble_stop_reason`
            !! above. Absent means no filtering (every Stop Condition allowed) -- deliberately
            !! nullable, not annotated with a generated default: the generator only evaluates
            !! constant *scalar* expressions for that annotation (`codegen_guide.md` section 5.5)
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible
            !! Per-ensemble eligibility from this criterion alone
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)

        call filter_ensembles_by_stop_condition(&
            n_ensembles = n_ensembles,&
            ensemble_stop_reason = ensemble_stop_reason,&
            allowed_stop_reasons = allowed_stop_reasons,&
            eligible = eligible,&
            ierr = ierr&
        )
    end subroutine filter_ensembles_by_stop_condition_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_filter(module):filter_ensembles_by_dimension(subroutine)]]
    !| `filter_dim_min <= d <= filter_dim_max`, both inclusive, each independently optional (an absent bound
    !| contributes no constraint on that side). An ensemble with no final accepted state at all
    !| (`ensemble_has_final(e)` false) is never eligible under this criterion once at least one
    !| of `filter_dim_min`/`filter_dim_max` is supplied -- there is no `d` to judge. Both bounds absent is a true
    !| no-op (every ensemble eligible, `ensemble_has_final` not even consulted), matching
    !| `filter_ensembles_by_stop_condition_impl`'s own "omitted means unconstrained"
    !| convention. `filter_dim_min` could in principle default to `0_int32` (a genuine constant
    !| expression), but is left nullable like `filter_dim_max` (whose own natural default, `n_dimensions`,
    !| is a runtime value and so cannot be a generated default) rather than have the two bounds
    !| of one range behave asymmetrically.
    subroutine filter_ensembles_by_dimension_c(&
            n_dimensions,&
            n_ensembles,&
            ensemble_d_final,&
            ensemble_has_final,&
            filter_dim_min,&
            filter_dim_max,&
            eligible,&
            ierr&
        ) bind(C, name="filter_ensembles_by_dimension_c")
        use tox_shape_truthful_clustering_filter, only: filter_ensembles_by_dimension

        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_d_final
            !! Each ensemble's final accepted intrinsic dimension, see
            !! `ensemble_final_observable`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        logical(c_bool), dimension(n_ensembles), intent(in), target :: ensemble_has_final
            !! Whether each ensemble has a final accepted state at all, see
            !! `ensemble_final_observable`
        integer(c_int), intent(in), optional :: filter_dim_min
            !! Minimum tolerated final intrinsic dimension, inclusive
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), optional :: filter_dim_max
            !! Maximum tolerated final intrinsic dimension, inclusive
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible
            !! Per-ensemble eligibility from this criterion alone
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_has_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)

        call filter_ensembles_by_dimension(&
            n_dimensions = n_dimensions,&
            n_ensembles = n_ensembles,&
            ensemble_d_final = ensemble_d_final,&
            ensemble_has_final = ensemble_has_final,&
            filter_dim_min = filter_dim_min,&
            filter_dim_max = filter_dim_max,&
            eligible = eligible,&
            ierr = ierr&
        )
    end subroutine filter_ensembles_by_dimension_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_filter(module):filter_ensembles_by_var_explained(subroutine)]]
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
    subroutine filter_ensembles_by_var_explained_c(&
            n_dimensions,&
            n_ensembles,&
            ensemble_S_final,&
            ensemble_d_final,&
            ensemble_k_final,&
            ensemble_has_final,&
            var_explained_min,&
            eligible,&
            ierr&
        ) bind(C, name="filter_ensembles_by_var_explained_c")
        use tox_shape_truthful_clustering_filter, only: filter_ensembles_by_var_explained

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_ensembles
            !! Number of ensembles N_E
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_dimensions, n_ensembles), intent(in), target :: ensemble_S_final
            !! Each ensemble's final accepted singular values, see `ensemble_final_observable`
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_d_final
            !! Each ensemble's final accepted intrinsic dimension
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_k_final
            !! Each ensemble's final accepted size
            !! The minimum valid value is `0_int32`.
        logical(c_bool), dimension(n_ensembles), intent(in), target :: ensemble_has_final
            !! Whether each ensemble has a final accepted state at all
        real(c_double), intent(in), optional :: var_explained_min
            !! Minimum tolerated fraction of variance explained by the tangent subspace,
            !! inclusive
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible
            !! Per-ensemble eligibility from this criterion alone
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_final, n_dimensions * n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(ensemble_has_final, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)

        call filter_ensembles_by_var_explained(&
            n_dimensions = n_dimensions,&
            n_ensembles = n_ensembles,&
            ensemble_S_final = ensemble_S_final,&
            ensemble_d_final = ensemble_d_final,&
            ensemble_k_final = ensemble_k_final,&
            ensemble_has_final = ensemble_has_final,&
            var_explained_min = var_explained_min,&
            eligible = eligible,&
            ierr = ierr&
        )
    end subroutine filter_ensembles_by_var_explained_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_filter(module):filter_ensembles(subroutine)]]
    !| Extracts each ensemble's final accepted state once (via `ensemble_final_observable`,
    !| shared so this and `tox_stc_json`'s own reporting never derive it two different ways),
    !| then calls each of this module's three per-criterion filters and combines their masks
    !| with a plain logical AND. Also returns the three individual masks, not just the
    !| combination -- so a caller (the report, in particular) can say *which* criterion excluded
    !| a given ensemble, not merely that one did. Supplying none of `allowed_stop_reasons`/
    !| `filter_dim_min`/`filter_dim_max`/`var_explained_min` makes every ensemble eligible (all four masks
    !| all-`.true.`), matching each individual filter's own no-op convention.
    subroutine filter_ensembles_c(&
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
            ensemble_stop_reason,&
            allowed_stop_reasons,&
            filter_dim_min,&
            filter_dim_max,&
            var_explained_min,&
            eligible,&
            eligible_by_stop_condition,&
            eligible_by_dimension,&
            eligible_by_var_explained,&
            ierr&
        ) bind(C, name="filter_ensembles_c")
        use tox_shape_truthful_clustering_filter, only: filter_ensembles

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
            !! merged output
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_S_history
            !! Per-ensemble trailing singular values
        real(c_double), dimension(n_dimensions, o, n_ensembles), intent(in), target :: ensemble_mu_history
            !! Per-ensemble trailing centers
        real(c_double), dimension(o, n_ensembles), intent(in), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        integer(c_int), dimension(o, n_ensembles), intent(in), target :: ensemble_k_history
            !! Per-ensemble trailing sizes
        logical(c_bool), dimension(o, n_ensembles), intent(in), target :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself accepted
        integer(c_int), dimension(n_ensembles), intent(in), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition, see `filter_ensembles_by_stop_condition_impl`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `4_int32`.
        logical(c_bool), dimension(4), intent(in), optional :: allowed_stop_reasons
            !! See `filter_ensembles_by_stop_condition_impl`
        integer(c_int), intent(in), optional :: filter_dim_min
            !! See `filter_ensembles_by_dimension_impl`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), optional :: filter_dim_max
            !! See `filter_ensembles_by_dimension_impl`
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), intent(in), optional :: var_explained_min
            !! See `filter_ensembles_by_var_explained_impl`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible
            !! Combined eligibility: `.true.` only where all three per-criterion masks are
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_stop_condition
            !! See `filter_ensembles_by_stop_condition_impl`
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_dimension
            !! See `filter_ensembles_by_dimension_impl`
        logical(c_bool), dimension(n_ensembles), intent(out), target :: eligible_by_var_explained
            !! See `filter_ensembles_by_var_explained_impl`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

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
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_stop_condition, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_dimension, n_ensembles)
        M_CHECK_ARRAY_NON_NULL(eligible_by_var_explained, n_ensembles)

        call filter_ensembles(&
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
            ensemble_stop_reason = ensemble_stop_reason,&
            allowed_stop_reasons = allowed_stop_reasons,&
            filter_dim_min = filter_dim_min,&
            filter_dim_max = filter_dim_max,&
            var_explained_min = var_explained_min,&
            eligible = eligible,&
            eligible_by_stop_condition = eligible_by_stop_condition,&
            eligible_by_dimension = eligible_by_dimension,&
            eligible_by_var_explained = eligible_by_var_explained,&
            ierr = ierr&
        )
    end subroutine filter_ensembles_c

end module tox_shape_truthful_clustering_filter_c
#endif
