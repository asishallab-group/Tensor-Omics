! Macros (M_ prefix, as in src/f42/serde/macros.h) for values that more than one module of the
! tox_data_integration family needs. A dummy extent or a DM_DEFAULT cannot name a Fortran
! parameter -- the generated C and R bindings need the number itself -- so each such value is a
! macro, written once, here, and every module that needs the literal includes this file.
! M_MAX_N_BINS, M_MAX_CANDIDATE_PAIRS and the M_MODE_BIN_COUNT_*, M_MODE_PLATEAU_*,
! M_MODE_NEIGHBORHOOD_SIZE_*, M_MODE_DIVERGENCE_* and M_METHOD_JOIN_* values also back a parameter
! of the module that owns them; the *_DEFAULT values have no parameter and are used directly. A
! value only one module uses stays in that module as a CM_ macro, except that a mode's values are
! kept together: once one of them is needed here, the whole mode family is here.

! Fixed ceiling on any histogram bin count: MAX_N_BINS of tox_data_integration_js_comp_test_binning_impl,
! and the extent of the bin-sized work arrays of the bin-count, search and final-run modules.
#define M_MAX_N_BINS 256

! `bin_count_mode`'s values (the MODE_BIN_COUNT_* of tox_data_integration_js_comp_test_binning_impl):
! how a reference point's histogram bin count is chosen. Stone's cross-validation with every study's
! score counted equally (M_MODE_BIN_COUNT_STONE_EQUAL) is the default of every routine that takes the
! mode; determine_bin_count_stone_impl rejects occupancy.
#define M_MODE_BIN_COUNT_OCCUPANCY 0_int32
#define M_MODE_BIN_COUNT_STONE_EQUAL 1_int32
#define M_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED 2_int32
#define M_MODE_BIN_COUNT_STONE_POOLED 3_int32

! Defaults of the occupancy bin-count search (determine_bin_count_occupancy_impl, also forwarded by
! the search and final-run modules and the study-label Permutation-Test). M_OCCUPANCY_GAMMA_DEFAULT (growth of a candidate bin count)
! is unrelated to GAMMA (decay of the fixed-k grid's n_points).
#define M_OCCUPANCY_M_MIN_DEFAULT 3_int32
#define M_OCCUPANCY_M_MAX_DEFAULT 120_int32
#define M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT 10_int32
#define M_OCCUPANCY_GAMMA_DEFAULT 1.25_real64
#define M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.05_real64
#define M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.95_real64

! Defaults of construct_adaptive_neighborhoods_impl's optional knobs (tox_data_integration_preprocessing_impl),
! shared by DM_DEFAULT and M_DEFAULT_VAL there, in the adaptive parameter search and in the adaptive Permutation-Test.
#define M_ADAPTIVE_TAU_DEFAULT 0.1_real64
#define M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT 1.0_real64
#define M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT 0_int32
#define M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT 1_int32

! Ceiling on the number of parameter-search candidates: MAX_CANDIDATE_PAIRS of
! tox_data_integration_js_comp_test_search_impl, and the extent of the candidate and trace arrays
! of the fixed-k and adaptive parameter searches.
#define M_MAX_CANDIDATE_PAIRS 16

! `neighborhood_size_mode`'s values (the MODE_NEIGHBORHOOD_SIZE_* of
! tox_data_integration_js_comp_test_fixed_k_search_impl): which neighbor-count factor of KX_FACTORS
! the fixed-k grid pairs with each n_points -- 0 the large neighborhoods (KX_FACTORS(1), the
! default), 1 the small ones (KX_FACTORS(2)). Here because the fixed-k search and the study-label
! Permutation-Test both name the default.
#define M_MODE_NEIGHBORHOOD_SIZE_LARGE 0_int32
#define M_MODE_NEIGHBORHOOD_SIZE_SMALL 1_int32

! Defaults of the fixed-k grid's two options (calc_fixed_k_n_points_range of
! tox_data_integration_js_comp_test_fixed_k_search_impl), shared by the DM_DEFAULT and
! M_DEFAULT_VAL of the candidate bounds, the grid, the fixed-k search and the study-label
! Permutation-Test: `grid_points_factor`, the factor in the grid's first n_points
! ceiling(factor*sqrt(G)), and `grid_min_points`, the grid's floor (the lower clamp of its first
! n_points and its last, coarsest candidate).
#define M_GRID_POINTS_FACTOR_DEFAULT 4.0_real64
#define M_GRID_MIN_POINTS_DEFAULT 200_int32

! `plateau_mode`'s values (the MODE_PLATEAU_* of tox_data_integration_js_comp_test_search_impl).
! The change of the JSCompTest values (effect size) is the DM_DEFAULT and M_DEFAULT_VAL of both
! parameter searches; CI overlap stays available as an option.
#define M_MODE_PLATEAU_CI_OVERLAP 0_int32
#define M_MODE_PLATEAU_EFFECT_SIZE 1_int32
#define M_MODE_PLATEAU_BOTH 2_int32

! `divergence_mode`'s values (the MODE_DIVERGENCE_* of tox_data_integration_jsd_impl): how each
! study's per-point divergence is computed -- 0 the JSD of the study against the consensus of all
! studies, 1 (the default) for exactly two studies the pairwise JSD of the two studies, shared by
! both, and for any other number of studies the same as 0. Here because the final runs, both
! parameter searches, the bootstrap and the study-label Permutation-Test all name the default.
#define M_MODE_DIVERGENCE_SYMMETRIC 0_int32
#define M_MODE_DIVERGENCE_PAIRWISE_K2 1_int32

! `join_method`'s values (the METHOD_JOIN_* of tox_data_integration_js_comp_test_search_impl): how
! the studies' confidence-interval overlaps are joined. Both parameter searches take it mandatory; the
! study-label Permutation-Test (tox_data_integration_js_comp_test_permutation_fixed_k_impl and
! tox_data_integration_js_comp_test_permutation_adaptive_impl) forwards it with M_METHOD_JOIN_MIN as its
! DM_DEFAULT and M_DEFAULT_VAL.
#define M_METHOD_JOIN_MIN 0_int32
#define M_METHOD_JOIN_MAX 1_int32
#define M_METHOD_JOIN_MEDIAN 2_int32

! Defaults of the two admissibility/plateau thresholds `min_neighbor_overlap` (gate 1) and
! `succeeding_ci_overlap`, shared by DM_DEFAULT and M_DEFAULT_VAL in both parameter searches and the Permutation-Test.
#define M_MIN_NEIGHBOR_OVERLAP_DEFAULT 0.1_real64
#define M_SUCCEEDING_CI_OVERLAP_DEFAULT 0.9_real64

! Defaults of the bootstrap (tox_data_integration_stats_impl), shared by its DM_DEFAULT and
! M_DEFAULT_VAL and by the DM_DEFAULT of every routine that forwards them: both parameter searches
! and the study-label Permutation-Test.
#define M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT 2.5_real64
#define M_RANDOM_SEED_DEFAULT 42_int32

! Default number of bootstraps of the study-label Permutation-Test's rerun parameter search
! (mandatory in both parameter searches themselves), shared by DM_DEFAULT and M_DEFAULT_VAL of the
! four Permutation-Test routines (the modules tox_data_integration_js_comp_test_permutation_fixed_k_impl
! and tox_data_integration_js_comp_test_permutation_adaptive_impl) and of their two sizing routines
! (tox_data_integration_js_comp_test_permutation_bounds_impl).
#define M_PERMUTATION_N_BOOTSTRAPS_DEFAULT 1000_int32

! Default of `rerun_parameter_search` of the four Permutation-Test routines: the parameter search is not rerun
! on every shuffled data set, shared by DM_DEFAULT and M_DEFAULT_VAL.
#define M_RERUN_PARAMETER_SEARCH_DEFAULT .false.

! Defaults of the relative-effect-size plateau criterion (provisional, to be validated
! empirically), shared by DM_DEFAULT and M_DEFAULT_VAL in both parameter searches and the Permutation-Test.
#define M_DELTA_MEDIAN_THRESHOLD_DEFAULT 0.05_real64
#define M_DELTA_MAX_THRESHOLD_DEFAULT 0.10_real64
#define M_DELTA_EPSILON_DEFAULT 1.0e-10_real64
#define M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT 2_int32
! Default of `delta_abs_tolerance` (check_effect_size_plateau_condition_impl, both parameter
! searches and the study-label Permutation-Test): changes of a study's JSD up to 0.0001 count as no
! change, so a JSD near 0 does not turn tiny absolute changes into large relative ones.
#define M_DELTA_ABS_TOLERANCE_DEFAULT 0.0001_real64

! The order condition of `lower_residual_range_quantile`, one doc sentence for every routine that
! takes both quantiles (the occupancy and Stone bin-count searches, both parameter searches,
! both final runs and the Permutation-Test), written once. It names no default: the docs of upper_residual_range_quantile
! show it.
#define M_LOWER_RESIDUAL_RANGE_QUANTILE_ORDER_DOC Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.

! The meaning of `delta_abs_tolerance`, one doc sentence for every routine that takes it
! (check_effect_size_plateau_condition_impl, both parameter searches and the four study-label
! Permutation-Test routines), written once. It names no default: the generated docs show it.
#define M_DELTA_ABS_TOLERANCE_DOC Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
