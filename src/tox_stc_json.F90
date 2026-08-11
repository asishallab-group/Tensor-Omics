#include <src/macros.h>

!> Serialization of Shape Truthful Clustering (STC) pipeline results into the JSON format
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
module tox_stc_json
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use f42_json, only: json_object, json_array, json_value, serialize_json_object
    use f42_utils, only: open_file
    use tox_errors, only: is_err, set_ok, validate_dimension_size, validate_in_range_int
    use tox_stc_html_assets, only: REPORT_TEMPLATE_HEAD, REPORT_TEMPLATE_MID, REPORT_TEMPLATE_TAIL, D3_JS
    use tox_shape_truthful_clustering_kernel, only: STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
        STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_REPORT, MODE_MERGE_OVERLAP_COEFFICIENT, &
        MODE_MERGE_ANY
    M_IMPLICIT_NONE

    private
    public :: serialize_stc_results_as_json
    public :: write_stc_interactive_html_report

    integer(int32), parameter :: STC_JSON_MAX_ENSEMBLE_KEYS = 13_int32
    integer(int32), parameter :: STC_JSON_MAX_POINT_KEYS = 7_int32
    integer(int32), parameter :: STC_JSON_MAX_PARAM_KEYS = 23_int32

contains

    !> M_EXPORT_C
    !| summary: Serializes an STC run's raw pipeline results as JSON
    !| AUTHOR_ASIS_HALLAB
    !| Matches the schema consumed by misc/STC-experiments/interactive_template.html's D3 report
    subroutine serialize_stc_results_as_json(filename, &
                                             n_dimensions, n_vectors, n_selected_seed, o, max_group_size, &
                                             n_super_ensembles, &
                                             vectors, dim_names, &
                                             seed_selection_mask, &
                                             ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
                                             ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                                             ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                                             ensemble_low_confidence_masks, &
                                             super_ensembles, &
                                             k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, &
                                             RMSE_change_max, f_max, a, exclusion_radius_percentile, &
                                             bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, &
                                             estimated_k_min, estimated_k_density, estimated_density_quantile, &
                                             estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, &
                                             estimated_d_max, &
                                             ierr)
        character(len=*), intent(in) :: filename
            !! Name of the JSON file to write
        integer(int32), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(int32), intent(in), target :: o
            !! Trailing observable-history window depth
        integer(int32), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(int32), intent(in) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        real(real64), intent(in), target :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        character(len=*), intent(in), target :: dim_names(n_dimensions)
            !! Per-dimension display name
        logical, intent(in) :: seed_selection_mask(n_vectors)
            !! Seed selection, see `seeds`
        logical, intent(in), target :: ensemble_masks(n_vectors, n_selected_seed)
            !! Per-ensemble accepted membership, one column per seed
        integer(int32), intent(in) :: ensemble_stop_reason(n_selected_seed)
            !! Per-ensemble Stop Condition
        real(real64), intent(in), target :: ensemble_growth_radii(n_selected_seed)
            !! Per-ensemble growth radius
        real(real64), intent(in), target :: ensemble_U_history(n_dimensions, n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing tangent+normal bases
        real(real64), intent(in), target :: ensemble_S_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing singular values
        integer(int32), intent(in), target :: ensemble_d_history(o, n_selected_seed)
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), intent(in), target :: ensemble_G_history(o, n_selected_seed)
            !! Per-ensemble trailing spectral gaps
        real(real64), intent(in), target :: ensemble_mu_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing centers
        integer(int32), intent(in) :: ensemble_k_history(o, n_selected_seed)
            !! Per-ensemble trailing sizes
        logical, intent(in) :: ensemble_low_confidence_masks(n_vectors, n_selected_seed)
            !! Per-ensemble iteration-1 fallback membership
        integer(int32), intent(in), target :: super_ensembles(max_group_size, n_selected_seed*(n_selected_seed - 1))
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(int32), intent(in), target :: k_min
            !! This run's neighborhood size for each seed's growth radius
        integer(int32), intent(in), target :: k_density
            !! This run's density estimation neighborhood size
        real(real64), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! This run's maximum tolerated chordal distance between tangent bases
        integer(int32), intent(in), target :: d_max
            !! This run's maximum tolerated change in intrinsic dimension
        real(real64), intent(in), target :: G_max
            !! This run's maximum tolerated |log(G_tp1/G_t)|
        real(real64), intent(in), target :: RMSE_change_max
            !! This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
        real(real64), intent(in), target :: f_max
            !! This run's ensemble size fraction of N above which growth is abandoned
        integer(int32), intent(in), target :: a
            !! This run's minimum accepted-iteration count for a stable rejection
        real(real64), intent(in), target :: exclusion_radius_percentile
            !! This run's seeding exclusion radius percentile
        real(real64), intent(in), target :: bandwidth_percentile
            !! This run's density-estimate kernel bandwidth percentile
        integer(int32), intent(in) :: reconciliation_mode
            !! This run's `ensemble_reconciliation` mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Report intersecting pairs only | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]] |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]] |
        real(real64), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
        integer(int32), intent(in), optional, target :: estimated_k_min
            !! `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
        integer(int32), intent(in), optional, target :: estimated_k_density
            !! `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
        real(real64), intent(in), optional, target :: estimated_density_quantile
            !! `estimate_stc_parameters`'s proposed density quantile, if estimation was used
        real(real64), intent(in), optional, target :: estimated_chordal_dist_max_as_prcnt_of_range
            !! `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
            !! estimation was used
        real(real64), intent(in), optional, target :: estimated_G_max
            !! `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
        integer(int32), intent(in), optional, target :: estimated_d_max
            !! `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: unit

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_seed, ierr, arg_pos=4_int32, min=0_int32)
        call validate_dimension_size(o, ierr, arg_pos=5_int32)
        call validate_in_range_int(max_group_size, ierr, arg_pos=6_int32, min=0_int32)
        call validate_in_range_int(n_super_ensembles, ierr, arg_pos=7_int32, min=0_int32)
        if (is_err(ierr)) return

        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return

        call stc_build_and_serialize_json(unit, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, &
            exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, &
            estimated_k_min, estimated_k_density, estimated_density_quantile, &
            estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, estimated_d_max)

        close (unit)
    end subroutine serialize_stc_results_as_json

    !> M_EXPORT_C
    !| summary: Writes a self-contained interactive HTML report for an STC run
    !| AUTHOR_ASIS_HALLAB
    !| Concatenates the vendored D3 library, the D3 report template (both baked in at compile
    !| time, see `tox_stc_html_assets`), and this run's own results as JSON into one file
    subroutine write_stc_interactive_html_report(filename, &
                                                 n_dimensions, n_vectors, n_selected_seed, o, max_group_size, &
                                                 n_super_ensembles, &
                                                 vectors, dim_names, &
                                                 seed_selection_mask, &
                                                 ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
                                                 ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                                                 ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                                                 ensemble_low_confidence_masks, &
                                                 super_ensembles, &
                                                 k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, &
                                                 RMSE_change_max, f_max, a, exclusion_radius_percentile, &
                                                 bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, &
                                                 estimated_k_min, estimated_k_density, estimated_density_quantile, &
                                                 estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, &
                                                 estimated_d_max, &
                                                 ierr)
        character(len=*), intent(in) :: filename
            !! Name of the HTML file to write
        integer(int32), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in), target :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(in), target :: n_selected_seed
            !! Number of selected seeds / accepted ensembles
        integer(int32), intent(in), target :: o
            !! Trailing observable-history window depth
        integer(int32), intent(in), target :: max_group_size
            !! Maximum number of ensembles one super-ensemble can hold
        integer(int32), intent(in) :: n_super_ensembles
            !! Number of leading columns of `super_ensembles` actually filled
        real(real64), intent(in), target :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        character(len=*), intent(in), target :: dim_names(n_dimensions)
            !! Per-dimension display name
        logical, intent(in) :: seed_selection_mask(n_vectors)
            !! Seed selection, see `seeds`
        logical, intent(in), target :: ensemble_masks(n_vectors, n_selected_seed)
            !! Per-ensemble accepted membership, one column per seed
        integer(int32), intent(in) :: ensemble_stop_reason(n_selected_seed)
            !! Per-ensemble Stop Condition
        real(real64), intent(in), target :: ensemble_growth_radii(n_selected_seed)
            !! Per-ensemble growth radius
        real(real64), intent(in), target :: ensemble_U_history(n_dimensions, n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing tangent+normal bases
        real(real64), intent(in), target :: ensemble_S_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing singular values
        integer(int32), intent(in), target :: ensemble_d_history(o, n_selected_seed)
            !! Per-ensemble trailing intrinsic dimensions
        real(real64), intent(in), target :: ensemble_G_history(o, n_selected_seed)
            !! Per-ensemble trailing spectral gaps
        real(real64), intent(in), target :: ensemble_mu_history(n_dimensions, o, n_selected_seed)
            !! Per-ensemble trailing centers
        integer(int32), intent(in) :: ensemble_k_history(o, n_selected_seed)
            !! Per-ensemble trailing sizes
        logical, intent(in) :: ensemble_low_confidence_masks(n_vectors, n_selected_seed)
            !! Per-ensemble iteration-1 fallback membership
        integer(int32), intent(in), target :: super_ensembles(max_group_size, n_selected_seed*(n_selected_seed - 1))
            !! One super-ensemble per column, 0-padded, see `ensemble_reconciliation`
        integer(int32), intent(in), target :: k_min
            !! This run's neighborhood size for each seed's growth radius
        integer(int32), intent(in), target :: k_density
            !! This run's density estimation neighborhood size
        real(real64), intent(in), target :: chordal_dist_max_as_prcnt_of_range
            !! This run's maximum tolerated chordal distance between tangent bases
        integer(int32), intent(in), target :: d_max
            !! This run's maximum tolerated change in intrinsic dimension
        real(real64), intent(in), target :: G_max
            !! This run's maximum tolerated |log(G_tp1/G_t)|
        real(real64), intent(in), target :: RMSE_change_max
            !! This run's maximum tolerated |log(RMSE_tp1/RMSE_t)|
        real(real64), intent(in), target :: f_max
            !! This run's ensemble size fraction of N above which growth is abandoned
        integer(int32), intent(in), target :: a
            !! This run's minimum accepted-iteration count for a stable rejection
        real(real64), intent(in), target :: exclusion_radius_percentile
            !! This run's seeding exclusion radius percentile
        real(real64), intent(in), target :: bandwidth_percentile
            !! This run's density-estimate kernel bandwidth percentile
        integer(int32), intent(in) :: reconciliation_mode
            !! This run's `ensemble_reconciliation` mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Report intersecting pairs only | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_REPORT(variable)]] |
            !! | Merge transitively at a minimum Overlap Coefficient | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_OVERLAP_COEFFICIENT(variable)]] |
            !! | Merge transitively on any intersection | [[tox_shape_truthful_clustering_reconciliation_kernel(module):MODE_MERGE_ANY(variable)]] |
        real(real64), intent(in), target :: min_overlap_coefficient
            !! This run's minimum Overlap Coefficient for `MODE_MERGE_OVERLAP_COEFFICIENT`
        integer(int32), intent(in), optional, target :: estimated_k_min
            !! `estimate_stc_parameters`'s proposed `k_min`, if estimation was used
        integer(int32), intent(in), optional, target :: estimated_k_density
            !! `estimate_stc_parameters`'s proposed `k_density`, if estimation was used
        real(real64), intent(in), optional, target :: estimated_density_quantile
            !! `estimate_stc_parameters`'s proposed density quantile, if estimation was used
        real(real64), intent(in), optional, target :: estimated_chordal_dist_max_as_prcnt_of_range
            !! `estimate_stc_parameters`'s proposed `chordal_dist_max_as_prcnt_of_range`, if
            !! estimation was used
        real(real64), intent(in), optional, target :: estimated_G_max
            !! `estimate_stc_parameters`'s proposed `G_max`, if estimation was used
        integer(int32), intent(in), optional, target :: estimated_d_max
            !! `estimate_stc_parameters`'s proposed `d_max`, if estimation was used
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        integer(int32) :: unit

        call set_ok(ierr)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_selected_seed, ierr, arg_pos=4_int32, min=0_int32)
        call validate_dimension_size(o, ierr, arg_pos=5_int32)
        call validate_in_range_int(max_group_size, ierr, arg_pos=6_int32, min=0_int32)
        call validate_in_range_int(n_super_ensembles, ierr, arg_pos=7_int32, min=0_int32)
        if (is_err(ierr)) return

        call open_file(filename, unit, .true., ierr)
        if (is_err(ierr)) return

        write (unit, "(A)", advance="no") REPORT_TEMPLATE_HEAD
        write (unit, "(A)", advance="no") D3_JS
        write (unit, "(A)", advance="no") REPORT_TEMPLATE_MID

        call stc_build_and_serialize_json(unit, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, RMSE_change_max, f_max, a, &
            exclusion_radius_percentile, bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, &
            estimated_k_min, estimated_k_density, estimated_density_quantile, &
            estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, estimated_d_max)

        write (unit, "(A)", advance="no") REPORT_TEMPLATE_TAIL

        close (unit)
    end subroutine write_stc_interactive_html_report

    !> Human-readable name for one of `ensemble_identification`'s `STOP_REASON_*` result codes.
    pure function stc_stop_reason_name(stop_reason) result(name)
        integer(int32), intent(in) :: stop_reason
        character(len=24) :: name

        select case (stop_reason)
            case (STOP_REASON_MAX_SIZE)
                name = 'max_size'
            case (STOP_REASON_REJECTED_AFTER_STABLE)
                name = 'rejected_after_stable'
            case (STOP_REASON_REJECTED_IMMEDIATELY)
                name = 'rejected_immediately'
            case (STOP_REASON_FIXED_POINT)
                name = 'fixed_point'
            case default
                name = 'error'
        end select
    end function stc_stop_reason_name

    !> Human-readable name for one of `ensemble_reconciliation`'s `MODE_*` values.
    pure function stc_reconciliation_mode_name(mode) result(name)
        integer(int32), intent(in) :: mode
        character(len=32) :: name

        select case (mode)
            case (MODE_REPORT)
                name = 'report'
            case (MODE_MERGE_OVERLAP_COEFFICIENT)
                name = 'merge_overlap_coefficient'
            case (MODE_MERGE_ANY)
                name = 'merge_any'
            case default
                name = 'unknown'
        end select
    end function stc_reconciliation_mode_name

    !> Column index of the last (most recently retained) populated iteration in one ensemble's
    !| trailing history, or 0 if the ensemble never produced one at all (only possible for
    !| `STOP_REASON_MAX_SIZE` firing at the bootstrap step itself, before any genuine SVD ever
    !| ran for that seed -- see `misc/mod_STC.md`, "Ensemble identification").
    pure function stc_last_history_index(k_history_col, o) result(idx)
        integer(int32), intent(in) :: o
        integer(int32), intent(in) :: k_history_col(o)
        integer(int32) :: idx

        integer(int32) :: j

        idx = 0
        do j = o, 1, -1
            if (k_history_col(j) /= 0) then
                idx = j
                exit
            end if
        end do
    end function stc_last_history_index

    !> Builds the full STC results document model and writes it, token by token, to `unit`.
    !| Not itself exported: both public entries call this with an already-open unit, differing
    !| only in what (if anything) surrounds the JSON in the file. Every array pointed at by the
    !| assembled tree is either one of this subroutine's own `target` dummy arguments or one of
    !| its local, automatic (explicit-shape) `target` workspace arrays -- both kinds live for
    !| the whole of this call, which is exactly as long as the tree itself is needed, since
    !| serialization happens before returning.
    subroutine stc_build_and_serialize_json(unit, &
                                            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, &
                                            n_super_ensembles, &
                                            vectors, dim_names, seed_selection_mask, &
                                            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
                                            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                                            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                                            ensemble_low_confidence_masks, super_ensembles, &
                                            k_min, k_density, chordal_dist_max_as_prcnt_of_range, d_max, G_max, &
                                            RMSE_change_max, f_max, a, exclusion_radius_percentile, &
                                            bandwidth_percentile, reconciliation_mode, min_overlap_coefficient, &
                                            estimated_k_min, estimated_k_density, estimated_density_quantile, &
                                            estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, &
                                            estimated_d_max)
        integer(int32), intent(in) :: unit
        integer(int32), intent(in), target :: n_dimensions
        integer(int32), intent(in), target :: n_vectors
        integer(int32), intent(in), target :: n_selected_seed
        integer(int32), intent(in), target :: o
        integer(int32), intent(in), target :: max_group_size
        integer(int32), intent(in) :: n_super_ensembles
        real(real64), intent(in), target :: vectors(n_dimensions, n_vectors)
        character(len=*), intent(in), target :: dim_names(n_dimensions)
        logical, intent(in) :: seed_selection_mask(n_vectors)
        logical, intent(in), target :: ensemble_masks(n_vectors, n_selected_seed)
        integer(int32), intent(in) :: ensemble_stop_reason(n_selected_seed)
        real(real64), intent(in), target :: ensemble_growth_radii(n_selected_seed)
        real(real64), intent(in), target :: ensemble_U_history(n_dimensions, n_dimensions, o, n_selected_seed)
        real(real64), intent(in), target :: ensemble_S_history(n_dimensions, o, n_selected_seed)
        integer(int32), intent(in), target :: ensemble_d_history(o, n_selected_seed)
        real(real64), intent(in), target :: ensemble_G_history(o, n_selected_seed)
        real(real64), intent(in), target :: ensemble_mu_history(n_dimensions, o, n_selected_seed)
        integer(int32), intent(in) :: ensemble_k_history(o, n_selected_seed)
        logical, intent(in) :: ensemble_low_confidence_masks(n_vectors, n_selected_seed)
        integer(int32), intent(in), target :: super_ensembles(max_group_size, n_selected_seed*(n_selected_seed - 1))
        integer(int32), intent(in), target :: k_min
        integer(int32), intent(in), target :: k_density
        real(real64), intent(in), target :: chordal_dist_max_as_prcnt_of_range
        integer(int32), intent(in), target :: d_max
        real(real64), intent(in), target :: G_max
        real(real64), intent(in), target :: RMSE_change_max
        real(real64), intent(in), target :: f_max
        integer(int32), intent(in), target :: a
        real(real64), intent(in), target :: exclusion_radius_percentile
        real(real64), intent(in), target :: bandwidth_percentile
        integer(int32), intent(in) :: reconciliation_mode
        real(real64), intent(in), target :: min_overlap_coefficient
        integer(int32), intent(in), optional, target :: estimated_k_min
        integer(int32), intent(in), optional, target :: estimated_k_density
        real(real64), intent(in), optional, target :: estimated_density_quantile
        real(real64), intent(in), optional, target :: estimated_chordal_dist_max_as_prcnt_of_range
        real(real64), intent(in), optional, target :: estimated_G_max
        integer(int32), intent(in), optional, target :: estimated_d_max

        ! -- params --------------------------------------------------------------------------
        type(json_object), target :: params_obj
        character(len=48), target :: param_keys(STC_JSON_MAX_PARAM_KEYS)
        type(json_value), target :: param_values(STC_JSON_MAX_PARAM_KEYS)
        character(len=32), target :: reconciliation_mode_name_buf
        integer(int32) :: n_params

        ! -- points ----------------------------------------------------------------------------
        type(json_array), target :: points_arr
        type(json_object), target :: point_objs(n_vectors)
        character(len=32), target :: point_keys(STC_JSON_MAX_POINT_KEYS)
        type(json_value), target :: point_values(STC_JSON_MAX_POINT_KEYS, n_vectors)
        integer(int32), target :: point_id_buf(n_vectors)
        integer(int32), target :: point_n_ensembles_buf(n_vectors)
        integer(int32), target :: point_n_low_conf_buf(n_vectors)
        integer(int32), target :: point_ensembles_buf(n_selected_seed, n_vectors)
        integer(int32) :: point_ensembles_count(n_vectors)
        integer(int32), target :: point_low_conf_buf(n_selected_seed, n_vectors)
        integer(int32) :: point_low_conf_count(n_vectors)
        integer(int32), target :: point_seed_of_buf(1, n_vectors)
        integer(int32) :: point_seed_of_count(n_vectors)
        type(json_array), target :: point_coords_arr(n_vectors)
        type(json_array), target :: point_ensembles_arr(n_vectors)
        type(json_array), target :: point_low_conf_arr(n_vectors)
        type(json_array), target :: point_seed_of_arr(n_vectors)

        ! -- ensembles ---------------------------------------------------------------------------
        type(json_array), target :: ensembles_arr
        type(json_object), target :: ensemble_objs(n_selected_seed)
        character(len=24), target :: ensemble_keys(STC_JSON_MAX_ENSEMBLE_KEYS, n_selected_seed)
        type(json_value), target :: ensemble_values(STC_JSON_MAX_ENSEMBLE_KEYS, n_selected_seed)
        integer(int32), target :: ensemble_id_buf(n_selected_seed)
        integer(int32), target :: ensemble_seed_point_id_buf(n_selected_seed)
        character(len=24), target :: ensemble_stop_reason_name_buf(n_selected_seed)
        integer(int32), target :: ensemble_size_buf(n_selected_seed)
        integer(int32), target :: ensemble_super_id_buf(n_selected_seed)
        type(json_array), target :: ensemble_mu_arr(n_selected_seed)
        type(json_array), target :: ensemble_u1_arr(n_selected_seed)
        type(json_array), target :: ensemble_u2_arr(n_selected_seed)

        ! -- super_ensembles ---------------------------------------------------------------------
        type(json_array), target :: super_ensembles_arr
        type(json_object), target :: se_objs(n_super_ensembles)
        character(len=16), target :: se_keys(2)
        type(json_value), target :: se_values(2, n_super_ensembles)
        integer(int32), target :: se_group_id_buf(n_super_ensembles)
        type(json_array), target :: se_ensemble_ids_arr(n_super_ensembles)

        ! -- overlap_coefficient_matrix ------------------------------------------------------
        type(json_array), target :: overlap_matrix_arr
        integer(int32), parameter :: OC_KEY_COUNT = 3_int32
        type(json_object), target :: oc_objs((n_selected_seed*(n_selected_seed - 1))/2)
        character(len=20), target :: oc_keys(OC_KEY_COUNT)
        type(json_value), target :: oc_values(OC_KEY_COUNT, (n_selected_seed*(n_selected_seed - 1))/2)
        integer(int32), target :: oc_a_buf((n_selected_seed*(n_selected_seed - 1))/2)
        integer(int32), target :: oc_b_buf((n_selected_seed*(n_selected_seed - 1))/2)
        real(real64), target :: oc_value_buf((n_selected_seed*(n_selected_seed - 1))/2)
        integer(int32) :: n_oc_pairs

        ! -- root ------------------------------------------------------------------------------
        type(json_object) :: root
        character(len=32), target :: root_keys(6)
        type(json_array), target :: root_dim_names
        type(json_value), target :: root_values(6)

        integer(int32) :: i, j, s, r, slot, idx, d_val, group_size, intersect_count
        integer(int32) :: vector_to_seed_index(n_vectors)

        ! ==== reconstruct the seed-index -> vector-index mapping =============================
        ! Mirrors ensemble_identification_merged_kernel's own construction of `seed_indices`
        ! exactly: column s of every per-ensemble output corresponds to the s-th .true. entry
        ! of seed_selection_mask, scanning vectors in order.
        vector_to_seed_index = 0
        s = 0
        do i = 1, n_vectors
            if (seed_selection_mask(i)) then
                s = s + 1
                ensemble_seed_point_id_buf(s) = i
                vector_to_seed_index(i) = s
            end if
        end do

        ! ==== params ===========================================================================
        n_params = 0
        n_params = n_params + 1; param_keys(n_params) = 'k_min'
        param_values(n_params)%value => k_min
        n_params = n_params + 1; param_keys(n_params) = 'k_density'
        param_values(n_params)%value => k_density
        n_params = n_params + 1; param_keys(n_params) = 'chordal_dist_max_as_prcnt_of_range'
        param_values(n_params)%value => chordal_dist_max_as_prcnt_of_range
        n_params = n_params + 1; param_keys(n_params) = 'd_max'
        param_values(n_params)%value => d_max
        n_params = n_params + 1; param_keys(n_params) = 'g_max'
        param_values(n_params)%value => G_max
        n_params = n_params + 1; param_keys(n_params) = 'rmse_change_max'
        param_values(n_params)%value => RMSE_change_max
        n_params = n_params + 1; param_keys(n_params) = 'f_max'
        param_values(n_params)%value => f_max
        n_params = n_params + 1; param_keys(n_params) = 'a'
        param_values(n_params)%value => a
        n_params = n_params + 1; param_keys(n_params) = 'o'
        param_values(n_params)%value => o
        n_params = n_params + 1; param_keys(n_params) = 'exclusion_radius_percentile'
        param_values(n_params)%value => exclusion_radius_percentile
        n_params = n_params + 1; param_keys(n_params) = 'bandwidth_percentile'
        param_values(n_params)%value => bandwidth_percentile
        n_params = n_params + 1; param_keys(n_params) = 'reconciliation_mode'
        reconciliation_mode_name_buf = stc_reconciliation_mode_name(reconciliation_mode)
        param_values(n_params)%value => reconciliation_mode_name_buf
        n_params = n_params + 1; param_keys(n_params) = 'min_overlap_coefficient'
        param_values(n_params)%value => min_overlap_coefficient
        n_params = n_params + 1; param_keys(n_params) = 'max_group_size'
        param_values(n_params)%value => max_group_size
        n_params = n_params + 1; param_keys(n_params) = 'n_vectors'
        param_values(n_params)%value => n_vectors
        n_params = n_params + 1; param_keys(n_params) = 'n_dimensions'
        param_values(n_params)%value => n_dimensions
        n_params = n_params + 1; param_keys(n_params) = 'n_ensembles'
        param_values(n_params)%value => n_selected_seed

        if (present(estimated_k_min)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_k_min'
            param_values(n_params)%value => estimated_k_min
        end if
        if (present(estimated_k_density)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_k_density'
            param_values(n_params)%value => estimated_k_density
        end if
        if (present(estimated_density_quantile)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_density_quantile'
            param_values(n_params)%value => estimated_density_quantile
        end if
        if (present(estimated_chordal_dist_max_as_prcnt_of_range)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_chordal_dist_max_as_prcnt_of_range'
            param_values(n_params)%value => estimated_chordal_dist_max_as_prcnt_of_range
        end if
        if (present(estimated_G_max)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_g_max'
            param_values(n_params)%value => estimated_G_max
        end if
        if (present(estimated_d_max)) then
            n_params = n_params + 1; param_keys(n_params) = 'estimated_d_max'
            param_values(n_params)%value => estimated_d_max
        end if

        params_obj%keys => param_keys(1:n_params)
        params_obj%values => param_values(1:n_params)

        ! ==== super_ensembles: derive each ensemble's group id first =========================
        ensemble_super_id_buf = 0
        do j = 1, n_super_ensembles
            do r = 1, max_group_size
                if (super_ensembles(r, j) == 0) cycle
                ensemble_super_id_buf(super_ensembles(r, j)) = j
            end do
        end do

        se_keys = [character(len=16) :: 'group_id', 'ensemble_ids']
        do j = 1, n_super_ensembles
            se_group_id_buf(j) = j
            group_size = count(super_ensembles(:, j) /= 0)
            if (group_size > 0) se_ensemble_ids_arr(j)%elements => super_ensembles(1:group_size, j)
            se_values(1, j)%value => se_group_id_buf(j)
            se_values(2, j)%value => se_ensemble_ids_arr(j)
            se_objs(j)%keys => se_keys
            se_objs(j)%values => se_values(:, j)
        end do
        if (n_super_ensembles > 0) super_ensembles_arr%elements => se_objs

        ! ==== ensembles ========================================================================
        do s = 1, n_selected_seed
            slot = 0

            slot = slot + 1; ensemble_keys(slot, s) = 'id'
            ensemble_id_buf(s) = s
            ensemble_values(slot, s)%value => ensemble_id_buf(s)

            slot = slot + 1; ensemble_keys(slot, s) = 'seed_point_id'
            ensemble_values(slot, s)%value => ensemble_seed_point_id_buf(s)

            slot = slot + 1; ensemble_keys(slot, s) = 'stop_reason'
            ensemble_stop_reason_name_buf(s) = stc_stop_reason_name(ensemble_stop_reason(s))
            ensemble_values(slot, s)%value => ensemble_stop_reason_name_buf(s)

            slot = slot + 1; ensemble_keys(slot, s) = 'growth_radius'
            ensemble_values(slot, s)%value => ensemble_growth_radii(s)

            slot = slot + 1; ensemble_keys(slot, s) = 'size'
            ensemble_size_buf(s) = count(ensemble_masks(:, s))
            ensemble_values(slot, s)%value => ensemble_size_buf(s)

            idx = stc_last_history_index(ensemble_k_history(:, s), o)
            if (idx > 0) then
                slot = slot + 1; ensemble_keys(slot, s) = 'd'
                ensemble_values(slot, s)%value => ensemble_d_history(idx, s)

                slot = slot + 1; ensemble_keys(slot, s) = 'G'
                ensemble_values(slot, s)%value => ensemble_G_history(idx, s)

                slot = slot + 1; ensemble_keys(slot, s) = 'mu'
                ensemble_mu_arr(s)%elements => ensemble_mu_history(:, idx, s)
                ensemble_values(slot, s)%value => ensemble_mu_arr(s)

                d_val = ensemble_d_history(idx, s)
                if (d_val >= 1) then
                    slot = slot + 1; ensemble_keys(slot, s) = 'u1'
                    ensemble_u1_arr(s)%elements => ensemble_U_history(:, 1, idx, s)
                    ensemble_values(slot, s)%value => ensemble_u1_arr(s)

                    slot = slot + 1; ensemble_keys(slot, s) = 's1'
                    ensemble_values(slot, s)%value => ensemble_S_history(1, idx, s)
                end if
                if (d_val >= 2) then
                    slot = slot + 1; ensemble_keys(slot, s) = 'u2'
                    ensemble_u2_arr(s)%elements => ensemble_U_history(:, 2, idx, s)
                    ensemble_values(slot, s)%value => ensemble_u2_arr(s)

                    slot = slot + 1; ensemble_keys(slot, s) = 's2'
                    ensemble_values(slot, s)%value => ensemble_S_history(2, idx, s)
                end if
            end if

            slot = slot + 1; ensemble_keys(slot, s) = 'super_ensemble_id'
            if (ensemble_super_id_buf(s) > 0) ensemble_values(slot, s)%value => ensemble_super_id_buf(s)

            ensemble_objs(s)%keys => ensemble_keys(1:slot, s)
            ensemble_objs(s)%values => ensemble_values(1:slot, s)
        end do
        if (n_selected_seed > 0) ensembles_arr%elements => ensemble_objs

        ! ==== points ===========================================================================
        point_keys = [character(len=32) :: 'id', 'coords', 'n_ensembles', 'n_low_confidence_ensembles', &
                     'ensembles', 'low_confidence_ensembles', 'seed_of']
        do i = 1, n_vectors
            point_id_buf(i) = i
            point_ensembles_count(i) = 0
            point_low_conf_count(i) = 0
            do s = 1, n_selected_seed
                if (ensemble_masks(i, s)) then
                    point_ensembles_count(i) = point_ensembles_count(i) + 1
                    point_ensembles_buf(point_ensembles_count(i), i) = s
                end if
                if (ensemble_low_confidence_masks(i, s)) then
                    point_low_conf_count(i) = point_low_conf_count(i) + 1
                    point_low_conf_buf(point_low_conf_count(i), i) = s
                end if
            end do
            point_n_ensembles_buf(i) = point_ensembles_count(i)
            point_n_low_conf_buf(i) = point_low_conf_count(i)

            point_seed_of_count(i) = 0
            if (seed_selection_mask(i)) then
                point_seed_of_count(i) = 1
                point_seed_of_buf(1, i) = vector_to_seed_index(i)
            end if

            point_coords_arr(i)%elements => vectors(:, i)
            if (point_ensembles_count(i) > 0) &
                point_ensembles_arr(i)%elements => point_ensembles_buf(1:point_ensembles_count(i), i)
            if (point_low_conf_count(i) > 0) &
                point_low_conf_arr(i)%elements => point_low_conf_buf(1:point_low_conf_count(i), i)
            if (point_seed_of_count(i) > 0) &
                point_seed_of_arr(i)%elements => point_seed_of_buf(1:1, i)

            point_values(1, i)%value => point_id_buf(i)
            point_values(2, i)%value => point_coords_arr(i)
            point_values(3, i)%value => point_n_ensembles_buf(i)
            point_values(4, i)%value => point_n_low_conf_buf(i)
            point_values(5, i)%value => point_ensembles_arr(i)
            point_values(6, i)%value => point_low_conf_arr(i)
            point_values(7, i)%value => point_seed_of_arr(i)

            point_objs(i)%keys => point_keys
            point_objs(i)%values => point_values(:, i)
        end do
        if (n_vectors > 0) points_arr%elements => point_objs

        ! ==== overlap_coefficient_matrix =======================================================
        oc_keys = [character(len=20) :: 'a', 'b', 'overlap_coefficient']
        n_oc_pairs = 0
        do i = 1, n_selected_seed - 1
            do j = i + 1, n_selected_seed
                intersect_count = count(ensemble_masks(:, i) .and. ensemble_masks(:, j))
                if (intersect_count < 1) cycle

                n_oc_pairs = n_oc_pairs + 1
                oc_a_buf(n_oc_pairs) = i
                oc_b_buf(n_oc_pairs) = j
                oc_value_buf(n_oc_pairs) = real(intersect_count, real64) / &
                    real(min(count(ensemble_masks(:, i)), count(ensemble_masks(:, j))), real64)

                oc_values(1, n_oc_pairs)%value => oc_a_buf(n_oc_pairs)
                oc_values(2, n_oc_pairs)%value => oc_b_buf(n_oc_pairs)
                oc_values(3, n_oc_pairs)%value => oc_value_buf(n_oc_pairs)

                oc_objs(n_oc_pairs)%keys => oc_keys
                oc_objs(n_oc_pairs)%values => oc_values(:, n_oc_pairs)
            end do
        end do
        if (n_oc_pairs > 0) overlap_matrix_arr%elements => oc_objs(1:n_oc_pairs)

        ! ==== root =============================================================================
        root_dim_names%elements => dim_names

        root_keys = [character(len=32) :: 'dim_names', 'params', 'points', 'ensembles', &
                    'super_ensembles', 'overlap_coefficient_matrix']
        root_values(1)%value => root_dim_names
        root_values(2)%value => params_obj
        root_values(3)%value => points_arr
        root_values(4)%value => ensembles_arr
        root_values(5)%value => super_ensembles_arr
        root_values(6)%value => overlap_matrix_arr
        root%keys => root_keys
        root%values => root_values

        call serialize_json_object(root, unit)
    end subroutine stc_build_and_serialize_json

end module tox_stc_json
