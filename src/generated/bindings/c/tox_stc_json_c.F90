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
!| `tox_shape_truthful_clustering_kernel`/`tox_shape_truthful_clustering_reconciliation_kernel`)
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
!| full pairwise Overlap Coefficient matrix) are computed here from the raw membership masks
!| -- `ensemble_reconciliation_kernel` itself only ever reports Overlap Coefficient along a
!| super-ensemble's own merge chain (`super_ensembles_overlap_coefficient`), never the full
!| N x N matrix the heatmap needs, so that matrix is recomputed directly from `ensemble_masks`
!| with the same `|intersect| / min(|A|,|B|)` formula, once per pair with a nonempty
!| intersection.
module tox_stc_json_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
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
            ensemble_low_confidence_masks,&
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
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="serialize_stc_results_as_json_c")
        use tox_stc_json, only: serialize_stc_results_as_json
        use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

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
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership
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
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
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
        character(len=:), allocatable :: filename_f
        character(len=:), allocatable, dimension(:) :: dim_names_f
        logical, dimension(n_vectors) :: seed_selection_mask_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_masks_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_low_confidence_masks_f
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
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_selected_seed*(n_selected_seed-1)))
        M_CHECK_ARRAY_NON_NULL(reconciliation_mode, 25)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(dim_names, dim_names_f, ierr)
        if (is_err(ierr)) return
        seed_selection_mask_f = seed_selection_mask
        ensemble_masks_f = ensemble_masks
        ensemble_low_confidence_masks_f = ensemble_low_confidence_masks
        block
            character(len=:), allocatable :: reconciliation_mode_f
            call c_char_1d_as_string(reconciliation_mode, reconciliation_mode_f, ierr)
            if (is_err(ierr)) return

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
            seed_selection_mask = seed_selection_mask_f,&
            ensemble_masks = ensemble_masks_f,&
            ensemble_stop_reason = ensemble_stop_reason,&
            ensemble_growth_radii = ensemble_growth_radii,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks_f,&
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
            ensemble_low_confidence_masks,&
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
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="write_stc_interactive_html_report_c")
        use tox_stc_json, only: write_stc_interactive_html_report
        use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_ANY, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_REPORT

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
        logical(c_bool), dimension(n_vectors, n_selected_seed), intent(in), target :: ensemble_low_confidence_masks
            !! Per-ensemble iteration-1 fallback membership
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
            !! | Mode                                                | Value                                                                                                    |
            !! |-----------------------------------------------------|----------------------------------------------------------------------------------------------------------|
            !! | Report intersecting pairs only                      | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]]                    |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection              | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]]                 |
        real(c_double), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
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
        character(len=:), allocatable :: filename_f
        character(len=:), allocatable, dimension(:) :: dim_names_f
        logical, dimension(n_vectors) :: seed_selection_mask_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_masks_f
        logical, dimension(n_vectors, n_selected_seed) :: ensemble_low_confidence_masks_f
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
        M_CHECK_ARRAY_NON_NULL(ensemble_low_confidence_masks, n_vectors * n_selected_seed)
        M_CHECK_ARRAY_NON_NULL(super_ensembles, max_group_size * (n_selected_seed*(n_selected_seed-1)))
        M_CHECK_ARRAY_NON_NULL(reconciliation_mode, 25)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call c_char_2d_as_string(dim_names, dim_names_f, ierr)
        if (is_err(ierr)) return
        seed_selection_mask_f = seed_selection_mask
        ensemble_masks_f = ensemble_masks
        ensemble_low_confidence_masks_f = ensemble_low_confidence_masks
        block
            character(len=:), allocatable :: reconciliation_mode_f
            call c_char_1d_as_string(reconciliation_mode, reconciliation_mode_f, ierr)
            if (is_err(ierr)) return

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
            seed_selection_mask = seed_selection_mask_f,&
            ensemble_masks = ensemble_masks_f,&
            ensemble_stop_reason = ensemble_stop_reason,&
            ensemble_growth_radii = ensemble_growth_radii,&
            ensemble_U_history = ensemble_U_history,&
            ensemble_S_history = ensemble_S_history,&
            ensemble_d_history = ensemble_d_history,&
            ensemble_G_history = ensemble_G_history,&
            ensemble_mu_history = ensemble_mu_history,&
            ensemble_k_history = ensemble_k_history,&
            ensemble_low_confidence_masks = ensemble_low_confidence_masks_f,&
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
