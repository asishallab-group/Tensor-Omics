#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_stc_json(module)]]
!| Serialization of Shape Truthful Clustering (STC) pipeline results into the JSON format
!| consumed by `misc/STC-experiments/interactive_template.html`'s D3 report, and into a
!| self-contained interactive HTML report combining that JSON with the report template and
!| the vendored D3 library (both baked in at compile time, see
!| `tox_stc_html_assets`/`helper/embed_stc_html_assets.py`). Builds a [[f42_json(module)]]
!| document model directly from STC's own raw result arrays (as returned by
!| `ensemble_identification_merged`/`ensemble_reconciliation`, see
!| `tox_shape_truthful_clustering_impl`/`tox_shape_truthful_clustering_reconciliation_impl`)
!| and writes it out -- this module is the STC-domain boundary on top of the generic
!| [[f42_json(module)]] serializer, the same role `tox_flyer_json` plays for the tox_flyer
!| viewer, and follows that module's exact single-subroutine pattern for the same reason: a
!| `json_object`/`json_array`'s pointer components must never outlive the frame that owns
!| their target storage, so tree-building and serialization happen inside one call.
!|
!| Deliberately writes with serial (`advance='no'`) writes directly into the destination
!| file's own stream unit, never materializing the JSON (or the assembled HTML) as an
!| in-memory string. One harmless, well-understood side effect: `close()` on a formatted
!| stream unit appends exactly one trailing newline beyond whatever was explicitly written --
!| irrelevant to a browser rendering HTML or a JSON parser, both of which ignore trailing
!| whitespace.
!|
!| `points`/`ensembles`/`super_ensembles` JSON keys not derivable from a single argument
!| (`id`, `n_ensembles`, per-point/per-ensemble membership lists, `super_ensemble_id`, the
!| full pairwise Overlap Coefficient matrix, per-point residual lengths, tangent line
!| endpoints, the two report-layer drift statistics below) are computed here from the raw
!| membership masks and history arrays -- `ensemble_reconciliation_impl` itself only ever
!| reports Overlap Coefficient along a super-ensemble's own merge chain
!| (`super_ensembles_overlap_coefficient`), never the full N x N matrix the heatmap needs, so
!| that matrix is recomputed directly from `ensemble_masks` with the same
!| `|intersect| / min(|A|,|B|)` formula, once per pair with a nonempty intersection.
!|
!| Two derived statistics reuse `tox_shape_truthful_clustering_accept_impl`'s own
!| `stc_chordal_distance` helper directly (made `public` there for exactly this reuse) rather
!| than re-deriving the formula: the **consecutive tangent-space drift** between each pair of
!| adjacent retained history columns (a genuine per-iteration quantity, fully reconstructable
!| from `ensemble_U_history` alone, unlike anything `accept_ensemble` itself tested), and the
!| **final accept-tested chordal distance** -- the one historical `accept_ensemble` criterion
!| (1) value that *is* exactly reconstructable after the fact, since the currently-stored
!| window minus its own last column, plus `ensemble_U_first`, is precisely the reference set
!| that was used to test the ensemble's actual final growth step. See `misc/mod_STC.md`,
!| "Ensemble Observable Plots", for the full rationale; neither statistic is stored as part of
!| `ensemble_identification`'s own output, keeping that SKG's kernels fully general and
!| iteration-unaware.
module tox_stc_json_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_f_pointer, c_int
    use, intrinsic :: iso_c_binding, only: c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: serialize_stc_results_as_json_c
    public :: write_stc_interactive_html_report_c

contains

    !> summary: C-wrapper for [[tox_stc_json(module):serialize_stc_results_as_json(subroutine)]]
    !| Matches the schema consumed by misc/STC-experiments/interactive_template.html's D3 report
    subroutine serialize_stc_results_as_json_c(&
            filename,&
            filename_strlen,&
            n_dimensions,&
            n_vectors,&
            n_selected_seed,&
            o,&
            max_group_size,&
            n_super_ensembles,&
            vectors,&
            dim_names,&
            dim_names_strlen,&
            seed_selection_mask,&
            ensemble_masks,&
            ensemble_stop_reason,&
            ensemble_growth_radii,&
            ensemble_U_history,&
            ensemble_S_history,&
            ensemble_d_history,&
            ensemble_G_history,&
            ensemble_mu_history,&
            ensemble_k_history,&
            ensemble_accepted_history,&
            ensemble_member_added_at_step,&
            ensemble_low_confidence_masks,&
            ensemble_U_first,&
            ensemble_d_first,&
            super_ensembles,&
            k_min,&
            k_density,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
            f_max,&
            a,&
            exclusion_radius_percentile,&
            bandwidth_percentile,&
            reconciliation_mode,&
            min_overlap_coefficient,&
            allowed_stop_reasons,&
            filter_dim_min,&
            filter_dim_max,&
            filter_var_explained_min,&
            ensemble_eligible,&
            ensemble_eligible_by_stop_condition,&
            ensemble_eligible_by_dimension,&
            ensemble_eligible_by_var_explained,&
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="serialize_stc_results_as_json_c")
        use tox_stc_json, only: serialize_stc_results_as_json
        use tox_shape_truthful_clustering_reconciliation_impl, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(c_int), intent(in), target :: dim_names_strlen
            !! length of the strings in `dim_names`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the JSON file to write
        integer(c_int), intent(in), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        character(len=1, kind=c_char), dimension(dim_names_strlen, n_dimensions), intent(in), target :: dim_names
            !! Per-dimension display name
        logical(c_bool), dimension(n_vectors), intent(in), target :: seed_selection_mask
            !! Seed selection, see `seeds`
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed
        integer(c_int), dimension(n_selected_seed), intent(in), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition
        real(c_double), dimension(n_selected_seed), intent(in), target :: ensemble_growth_radii
            !! Per-ensemble growth radius
        real(c_double), dimension(n_dimensions, n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_S_history
            !! Per-ensemble trailing singular values
        integer(c_int), dimension(o, n_selected_seed), intent(in), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(c_double), dimension(o, n_selected_seed), intent(in), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_mu_history
            !! Per-ensemble trailing centers
        integer(c_int), dimension(o, n_selected_seed), intent(in), target :: ensemble_k_history
            !! Per-ensemble trailing sizes
        logical(c_bool), dimension(o, n_selected_seed), intent(in), target :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself accepted
            !! -- `stc_push_ensemble_history` also pushes a *rejected* final candidate before
            !! `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
            !! `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is not always
            !! the ensemble's actual last accepted state; this module uses this array to find
            !! the last column that genuinely is (see `stc_last_accepted_history_index`)
        integer(c_int), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_member_added_at_step
            !! Per-ensemble growth-iteration-joined bookkeeping, see `ensemble_identification`'s
            !! `member_added_at_step`; this module only ever reads its column max (= T, the
            !! final accepted growth iteration), not the per-vector values themselves
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership
        real(c_double), dimension(n_dimensions, n_dimensions, n_selected_seed), intent(in), target :: ensemble_U_first
            !! Per-ensemble tangent+normal basis at the bootstrap iteration (iteration 1)
        integer(c_int), dimension(n_selected_seed), intent(in), target :: ensemble_d_first
            !! Per-ensemble intrinsic dimension at the bootstrap iteration
        integer(c_int), dimension(max_group_size, n_selected_seed*(n_selected_seed-1)), intent(in), target :: super_ensembles
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(c_int), intent(in), target :: k_min
            !! This run's neighborhood size for each seed's growth radius
        integer(c_int), intent(in), target :: k_density
            !! This run's density estimation neighborhood size
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! This run's maximum tolerated chordal distance between tangent bases
        integer(c_int), intent(in), target :: d_max
            !! This run's maximum tolerated change in intrinsic dimension
        real(c_double), intent(in), target :: G_max
            !! This run's maximum tolerated |log(G_tp1/G_t)|
        real(c_double), intent(in), target :: RMSE_change_max
            !! This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
        real(c_double), intent(in), target :: f_max
            !! This run's ensemble size fraction of N above which growth is abandoned
        integer(c_int), intent(in), target :: a
            !! This run's minimum accepted-iteration count for a stable rejection
        real(c_double), intent(in), target :: exclusion_radius_percentile
            !! This run's seeding exclusion radius percentile
        real(c_double), intent(in), target :: bandwidth_percentile
            !! This run's density-estimate kernel bandwidth percentile
        character(len=1, kind=c_char), dimension(25), intent(in), target :: reconciliation_mode
            !! This run's `ensemble_reconciliation` mode
            !!
            !! | Mode                                                | Value                                                                                                  |
            !! |-----------------------------------------------------|--------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_ANY(variable)]]                 |
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
        logical(c_bool), dimension(4), intent(in), optional :: allowed_stop_reasons
            !! This run's per-Stop-Condition eligibility actually used by
            !! `ensemble_reconciliation` -- reported here (as `params.excluded_stop_reasons`)
            !! for transparency only; this module no longer derives eligibility from it itself,
            !! see `ensemble_eligible` below
        integer(c_int), intent(in), optional :: filter_dim_min
            !! This run's minimum tolerated final intrinsic dimension for reconciliation
            !! eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `d_min`; reported for transparency only, same as `allowed_stop_reasons` above
        integer(c_int), intent(in), optional :: filter_dim_max
            !! This run's maximum tolerated final intrinsic dimension for reconciliation
            !! eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `d_max`; reported for transparency only, same as `allowed_stop_reasons` above
        real(c_double), intent(in), optional :: filter_var_explained_min
            !! This run's minimum tolerated final variance explained for reconciliation
            !! eligibility -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `var_explained_min`; reported for transparency only, same as
            !! `allowed_stop_reasons` above
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible
            !! Per-ensemble combined reconciliation eligibility actually used by
            !! `ensemble_reconciliation`, see its own `eligible` output
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_stop_condition
            !! See `ensemble_reconciliation`'s own `eligible_by_stop_condition`
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_dimension
            !! See `ensemble_reconciliation`'s own `eligible_by_dimension`
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_var_explained
            !! See `ensemble_reconciliation`'s own `eligible_by_var_explained`
        integer(c_int), intent(in), optional :: estimated_k_min
            !! `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
        integer(c_int), intent(in), optional :: estimated_k_density
            !! `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
        real(c_double), intent(in), optional :: estimated_density_quantile
            !! `estimate_stc_parameters`'s proposed density quantile, if estimation was used
        real(c_double), intent(in), optional :: estimated_chordal_dist_max_as_prcnt_of_range
            !! `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
            !! estimation was used
        real(c_double), intent(in), optional :: estimated_G_max
            !! `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
        integer(c_int), intent(in), optional :: estimated_d_max
            !! `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        character(len=filename_strlen), pointer :: filename_f
        character(len=dim_names_strlen), pointer, dimension(:) :: dim_names_f
        integer(int32) :: reconciliation_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_seed)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_NON_NULL(dim_names_strlen)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(f_max)
        M_CHECK_NON_NULL(a)
        M_CHECK_NON_NULL(exclusion_radius_percentile)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dim_names, dim_names_strlen * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(seed_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_growth_radii, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_history, n_dimensions * n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_accepted_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_member_added_at_step, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_first, n_dimensions * n_dimensions * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_first, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_selected_seed*(n_selected_seed-1)))
        M_CHECK_ARRAY_NON_NULL(reconciliation_mode, 25)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_stop_condition, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_dimension, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_var_explained, n_selected_seed)

        call c_f_pointer(c_loc(filename), filename_f)
        call c_f_pointer(c_loc(dim_names), dim_names_f, [n_dimensions])
        block
            character(len=:), pointer :: reconciliation_mode_f
            reconciliation_mode_f => c_char_as_view(reconciliation_mode)

            select case (reconciliation_mode_f)
                case ("report")
                    reconciliation_mode_mode_f = MODE_REPORT
                case ("merge_overlap_coefficient")
                    reconciliation_mode_mode_f = MODE_MERGE_OVERLAP_COEFFICIENT
                case ("merge_any")
                    reconciliation_mode_mode_f = MODE_MERGE_ANY
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call serialize_stc_results_as_json(&
            filename = filename_f,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            n_selected_seed = n_selected_seed,&
            o = o,&
            max_group_size = max_group_size,&
            n_super_ensembles = n_super_ensembles,&
            vectors = vectors,&
            dim_names = dim_names_f,&
            seed_selection_mask = seed_selection_mask,&
            ensemble_masks = ensemble_masks,&
            ensemble_stop_reason = ensemble_stop_reason,&
            ensemble_growth_radii = ensemble_growth_radii,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history,&
            ensemble_member_added_at_step = ensemble_member_added_at_step,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks,&
            ensemble_U_first = ensemble_U_first,&
            ensemble_d_first = ensemble_d_first,&
            super_ensembles = super_ensembles,&
            k_min = k_min,&
            k_density = k_density,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
            f_max = f_max,&
            a = a,&
            exclusion_radius_percentile = exclusion_radius_percentile,&
            bandwidth_percentile = bandwidth_percentile,&
            reconciliation_mode = reconciliation_mode_mode_f,&
            min_overlap_coefficient = min_overlap_coefficient,&
            allowed_stop_reasons = allowed_stop_reasons,&
            filter_dim_min = filter_dim_min,&
            filter_dim_max = filter_dim_max,&
            filter_var_explained_min = filter_var_explained_min,&
            ensemble_eligible = ensemble_eligible,&
            ensemble_eligible_by_stop_condition = ensemble_eligible_by_stop_condition,&
            ensemble_eligible_by_dimension = ensemble_eligible_by_dimension,&
            ensemble_eligible_by_var_explained = ensemble_eligible_by_var_explained,&
            estimated_k_min = estimated_k_min,&
            estimated_k_density = estimated_k_density,&
            estimated_density_quantile = estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range = estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max = estimated_G_max,&
            estimated_d_max = estimated_d_max,&
            ierr = ierr&
        )
    end subroutine serialize_stc_results_as_json_c

    !> summary: C-wrapper for [[tox_stc_json(module):write_stc_interactive_html_report(subroutine)]]
    !| Concatenates the vendored D3 library, the D3 report template (both baked in at compile
    !| time, see `tox_stc_html_assets`), and this run's own results as JSON into one file
    subroutine write_stc_interactive_html_report_c(&
            filename,&
            filename_strlen,&
            n_dimensions,&
            n_vectors,&
            n_selected_seed,&
            o,&
            max_group_size,&
            n_super_ensembles,&
            vectors,&
            dim_names,&
            dim_names_strlen,&
            seed_selection_mask,&
            ensemble_masks,&
            ensemble_stop_reason,&
            ensemble_growth_radii,&
            ensemble_U_history,&
            ensemble_S_history,&
            ensemble_d_history,&
            ensemble_G_history,&
            ensemble_mu_history,&
            ensemble_k_history,&
            ensemble_accepted_history,&
            ensemble_member_added_at_step,&
            ensemble_low_confidence_masks,&
            ensemble_U_first,&
            ensemble_d_first,&
            super_ensembles,&
            k_min,&
            k_density,&
            chordal_dist_max_as_prcnt_of_range,&
            d_max,&
            G_max,&
            RMSE_change_max,&
            f_max,&
            a,&
            exclusion_radius_percentile,&
            bandwidth_percentile,&
            reconciliation_mode,&
            min_overlap_coefficient,&
            allowed_stop_reasons,&
            filter_dim_min,&
            filter_dim_max,&
            filter_var_explained_min,&
            ensemble_eligible,&
            ensemble_eligible_by_stop_condition,&
            ensemble_eligible_by_dimension,&
            ensemble_eligible_by_var_explained,&
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="write_stc_interactive_html_report_c")
        use tox_stc_json, only: write_stc_interactive_html_report
        use tox_shape_truthful_clustering_reconciliation_impl, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(c_int), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(c_int), intent(in), target :: o
            !! Trailing observable-history window depth
        integer(c_int), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(c_int), intent(in), target :: dim_names_strlen
            !! length of the strings in `dim_names`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the HTML file to write
        integer(c_int), intent(in), target :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        character(len=1, kind=c_char), dimension(dim_names_strlen, n_dimensions), intent(in), target :: dim_names
            !! Per-dimension display name
        logical(c_bool), dimension(n_vectors), intent(in), target :: seed_selection_mask
            !! Seed selection, see `seeds`
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_masks
            !! Per-ensemble accepted membership, one column per seed
        integer(c_int), dimension(n_selected_seed), intent(in), target :: ensemble_stop_reason
            !! Per-ensemble Stop Condition
        real(c_double), dimension(n_selected_seed), intent(in), target :: ensemble_growth_radii
            !! Per-ensemble growth radius
        real(c_double), dimension(n_dimensions, n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_U_history
            !! Per-ensemble trailing tangent+normal bases
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_S_history
            !! Per-ensemble trailing singular values
        integer(c_int), dimension(o, n_selected_seed), intent(in), target :: ensemble_d_history
            !! Per-ensemble trailing intrinsic dimensions
        real(c_double), dimension(o, n_selected_seed), intent(in), target :: ensemble_G_history
            !! Per-ensemble trailing spectral gaps
        real(c_double), dimension(n_dimensions, o, n_selected_seed), intent(in), target :: ensemble_mu_history
            !! Per-ensemble trailing centers
        integer(c_int), dimension(o, n_selected_seed), intent(in), target :: ensemble_k_history
            !! Per-ensemble trailing sizes
        logical(c_bool), dimension(o, n_selected_seed), intent(in), target :: ensemble_accepted_history
            !! Whether the growth iteration retained in each history column was itself accepted
            !! -- `stc_push_ensemble_history` also pushes a *rejected* final candidate before
            !! `ensemble_identification` halts growth via `STOP_REASON_REJECTED_IMMEDIATELY`/
            !! `STOP_REASON_REJECTED_AFTER_STABLE`, so the last populated column is not always
            !! the ensemble's actual last accepted state; this module uses this array to find
            !! the last column that genuinely is (see `stc_last_accepted_history_index`)
        integer(c_int), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_member_added_at_step
            !! Per-ensemble growth-iteration-joined bookkeeping, see `ensemble_identification`'s
            !! `member_added_at_step`; this module only ever reads its column max (= T, the
            !! final accepted growth iteration), not the per-vector values themselves
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership
        real(c_double), dimension(n_dimensions, n_dimensions, n_selected_seed), intent(in), target :: ensemble_U_first
            !! Per-ensemble tangent+normal basis at the bootstrap iteration (iteration 1)
        integer(c_int), dimension(n_selected_seed), intent(in), target :: ensemble_d_first
            !! Per-ensemble intrinsic dimension at the bootstrap iteration
        integer(c_int), dimension(max_group_size, n_selected_seed*(n_selected_seed-1)), intent(in), target :: super_ensembles
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(c_int), intent(in), target :: k_min
            !! This run's neighborhood size for each seed's growth radius
        integer(c_int), intent(in), target :: k_density
            !! This run's density estimation neighborhood size
        real(c_double), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! This run's maximum tolerated chordal distance between tangent bases
        integer(c_int), intent(in), target :: d_max
            !! This run's maximum tolerated change in intrinsic dimension
        real(c_double), intent(in), target :: G_max
            !! This run's maximum tolerated |log(G_tp1/G_t)|
        real(c_double), intent(in), target :: RMSE_change_max
            !! This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
        real(c_double), intent(in), target :: f_max
            !! This run's ensemble size fraction of N above which growth is abandoned
        integer(c_int), intent(in), target :: a
            !! This run's minimum accepted-iteration count for a stable rejection
        real(c_double), intent(in), target :: exclusion_radius_percentile
            !! This run's seeding exclusion radius percentile
        real(c_double), intent(in), target :: bandwidth_percentile
            !! This run's density-estimate kernel bandwidth percentile
        character(len=1, kind=c_char), dimension(25), intent(in), target :: reconciliation_mode
            !! This run's `ensemble_reconciliation` mode
            !!
            !! | Mode                                                | Value                                                                                                  |
            !! |-----------------------------------------------------|--------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_impl(module):MODE_MERGE_ANY(variable)]]                 |
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
        logical(c_bool), dimension(4), intent(in), optional :: allowed_stop_reasons
            !! This run's per-Stop-Condition eligibility actually used by
            !! `ensemble_reconciliation` -- reported here (as `params.excluded_stop_reasons`)
            !! for transparency only; this module no longer derives eligibility from it itself,
            !! see `ensemble_eligible` below
        integer(c_int), intent(in), optional :: filter_dim_min
            !! This run's minimum tolerated final intrinsic dimension for reconciliation
            !! eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `d_min`; reported for transparency only, same as `allowed_stop_reasons` above
        integer(c_int), intent(in), optional :: filter_dim_max
            !! This run's maximum tolerated final intrinsic dimension for reconciliation
            !! eligibility, inclusive -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `d_max`; reported for transparency only, same as `allowed_stop_reasons` above
        real(c_double), intent(in), optional :: filter_var_explained_min
            !! This run's minimum tolerated final variance explained for reconciliation
            !! eligibility -- see `tox_shape_truthful_clustering_filter_impl`'s own
            !! `var_explained_min`; reported for transparency only, same as
            !! `allowed_stop_reasons` above
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible
            !! Per-ensemble combined reconciliation eligibility actually used by
            !! `ensemble_reconciliation`, see its own `eligible` output
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_stop_condition
            !! See `ensemble_reconciliation`'s own `eligible_by_stop_condition`
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_dimension
            !! See `ensemble_reconciliation`'s own `eligible_by_dimension`
        logical(c_bool), dimension(n_selected_seed), intent(in), target :: ensemble_eligible_by_var_explained
            !! See `ensemble_reconciliation`'s own `eligible_by_var_explained`
        integer(c_int), intent(in), optional :: estimated_k_min
            !! `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
        integer(c_int), intent(in), optional :: estimated_k_density
            !! `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
        real(c_double), intent(in), optional :: estimated_density_quantile
            !! `estimate_stc_parameters`'s proposed density quantile, if estimation was used
        real(c_double), intent(in), optional :: estimated_chordal_dist_max_as_prcnt_of_range
            !! `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
            !! estimation was used
        real(c_double), intent(in), optional :: estimated_G_max
            !! `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
        integer(c_int), intent(in), optional :: estimated_d_max
            !! `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        character(len=filename_strlen), pointer :: filename_f
        character(len=dim_names_strlen), pointer, dimension(:) :: dim_names_f
        integer(int32) :: reconciliation_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_seed)
        M_CHECK_NON_NULL(o)
        M_CHECK_NON_NULL(max_group_size)
        M_CHECK_NON_NULL(n_super_ensembles)
        M_CHECK_NON_NULL(dim_names_strlen)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(d_max)
        M_CHECK_NON_NULL(G_max)
        M_CHECK_NON_NULL(RMSE_change_max)
        M_CHECK_NON_NULL(f_max)
        M_CHECK_NON_NULL(a)
        M_CHECK_NON_NULL(exclusion_radius_percentile)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_NON_NULL(min_overlap_coefficient)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(dim_names, dim_names_strlen * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(seed_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(ensemble_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_stop_reason, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_growth_radii, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_history, n_dimensions * n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_S_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_G_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_mu_history, n_dimensions * o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_k_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_accepted_history, o * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_member_added_at_step, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_U_first, n_dimensions * n_dimensions * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_d_first, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_selected_seed*(n_selected_seed-1)))
        M_CHECK_ARRAY_NON_NULL(reconciliation_mode, 25)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_stop_condition, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_dimension, n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(ensemble_eligible_by_var_explained, n_selected_seed)

        call c_f_pointer(c_loc(filename), filename_f)
        call c_f_pointer(c_loc(dim_names), dim_names_f, [n_dimensions])
        block
            character(len=:), pointer :: reconciliation_mode_f
            reconciliation_mode_f => c_char_as_view(reconciliation_mode)

            select case (reconciliation_mode_f)
                case ("report")
                    reconciliation_mode_mode_f = MODE_REPORT
                case ("merge_overlap_coefficient")
                    reconciliation_mode_mode_f = MODE_MERGE_OVERLAP_COEFFICIENT
                case ("merge_any")
                    reconciliation_mode_mode_f = MODE_MERGE_ANY
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call write_stc_interactive_html_report(&
            filename = filename_f,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            n_selected_seed = n_selected_seed,&
            o = o,&
            max_group_size = max_group_size,&
            n_super_ensembles = n_super_ensembles,&
            vectors = vectors,&
            dim_names = dim_names_f,&
            seed_selection_mask = seed_selection_mask,&
            ensemble_masks = ensemble_masks,&
            ensemble_stop_reason = ensemble_stop_reason,&
            ensemble_growth_radii = ensemble_growth_radii,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_accepted_history = ensemble_accepted_history,&
            ensemble_member_added_at_step = ensemble_member_added_at_step,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks,&
            ensemble_U_first = ensemble_U_first,&
            ensemble_d_first = ensemble_d_first,&
            super_ensembles = super_ensembles,&
            k_min = k_min,&
            k_density = k_density,&
            chordal_dist_max_as_prcnt_of_range = chordal_dist_max_as_prcnt_of_range,&
            d_max = d_max,&
            G_max = G_max,&
            RMSE_change_max = RMSE_change_max,&
            f_max = f_max,&
            a = a,&
            exclusion_radius_percentile = exclusion_radius_percentile,&
            bandwidth_percentile = bandwidth_percentile,&
            reconciliation_mode = reconciliation_mode_mode_f,&
            min_overlap_coefficient = min_overlap_coefficient,&
            allowed_stop_reasons = allowed_stop_reasons,&
            filter_dim_min = filter_dim_min,&
            filter_dim_max = filter_dim_max,&
            filter_var_explained_min = filter_var_explained_min,&
            ensemble_eligible = ensemble_eligible,&
            ensemble_eligible_by_stop_condition = ensemble_eligible_by_stop_condition,&
            ensemble_eligible_by_dimension = ensemble_eligible_by_dimension,&
            ensemble_eligible_by_var_explained = ensemble_eligible_by_var_explained,&
            estimated_k_min = estimated_k_min,&
            estimated_k_density = estimated_k_density,&
            estimated_density_quantile = estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range = estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max = estimated_G_max,&
            estimated_d_max = estimated_d_max,&
            ierr = ierr&
        )
    end subroutine write_stc_interactive_html_report_c

end module tox_stc_json_c
#endif
