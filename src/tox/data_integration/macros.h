! Macros (M_ prefix, as in src/f42/serde/macros.h) for values that more than one module of the
! tox_data_integration family needs. A dummy extent or a DM_DEFAULT cannot name a Fortran
! parameter -- the generated C and R bindings need the number itself -- so each such value is a
! macro, written once, here, and every module that needs the literal includes this file.
! M_MAX_N_BINS, the M_MODE_BIN_COUNT_* values, M_MAX_CANDIDATE_PAIRS and M_MODE_PLATEAU_CI_OVERLAP
! also back a parameter of the module that owns them; the *_DEFAULT values have no parameter and
! are used directly. A value only one module uses stays in that module as a CM_ macro.

! Fixed ceiling on any histogram bin count: MAX_N_BINS of tox_data_integration_js_comp_test_binning_impl,
! and the extent of the bin-sized work arrays of the bin-count, search and final-run modules.
#define M_MAX_N_BINS 256

! `bin_count_mode`'s values (the MODE_BIN_COUNT_* of tox_data_integration_js_comp_test_binning_impl):
! how a reference point's histogram bin count is chosen. Occupancy is the default of all four
! pipeline entry points; determine_bin_count_stone_impl rejects it.
#define M_MODE_BIN_COUNT_OCCUPANCY 0_int32
#define M_MODE_BIN_COUNT_STONE_EQUAL 1_int32
#define M_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED 2_int32
#define M_MODE_BIN_COUNT_STONE_POOLED 3_int32

! Defaults of the occupancy bin-count search (determine_bin_count_occupancy_impl, also forwarded by
! the search and final-run modules). M_OCCUPANCY_GAMMA_DEFAULT (growth of a candidate bin count)
! is unrelated to GAMMA (decay of the fixed-k grid's n_points).
#define M_OCCUPANCY_M_MIN_DEFAULT 3_int32
#define M_OCCUPANCY_M_MAX_DEFAULT 120_int32
#define M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT 10_int32
#define M_OCCUPANCY_GAMMA_DEFAULT 1.25_real64
#define M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.05_real64
#define M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.95_real64

! Defaults of construct_adaptive_neighborhoods_impl's optional knobs (tox_data_integration_preprocessing_impl),
! shared by DM_DEFAULT and M_DEFAULT_VAL there and in the adaptive parameter search.
#define M_ADAPTIVE_TAU_DEFAULT 0.1_real64
#define M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT 1.0_real64
#define M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT 0_int32
#define M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT 1_int32

! Ceiling on the number of parameter-search candidates: MAX_CANDIDATE_PAIRS of
! tox_data_integration_js_comp_test_search_impl, and the extent of the candidate and trace arrays
! of the fixed-k and adaptive parameter searches.
#define M_MAX_CANDIDATE_PAIRS 16

! `plateau_mode`'s default, MODE_PLATEAU_CI_OVERLAP of tox_data_integration_js_comp_test_search_impl:
! the DM_DEFAULT and M_DEFAULT_VAL of both parameter searches.
#define M_MODE_PLATEAU_CI_OVERLAP 0_int32

! Defaults of the two admissibility/plateau thresholds `min_neighbor_overlap` (gate 1) and
! `succeeding_ci_overlap`, shared by DM_DEFAULT and M_DEFAULT_VAL in both parameter searches.
#define M_MIN_NEIGHBOR_OVERLAP_DEFAULT 0.1_real64
#define M_SUCCEEDING_CI_OVERLAP_DEFAULT 0.9_real64

! Defaults of the bootstrap and the permutation test (tox_data_integration_stats_impl), shared by
! their DM_DEFAULT and M_DEFAULT_VAL and by the DM_DEFAULT of every routine that forwards them:
! both parameter searches and both final runs.
#define M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT 2.5_real64
#define M_RANDOM_SEED_DEFAULT 42_int32

! Defaults of the relative-effect-size plateau criterion (provisional, to be validated
! empirically), shared by DM_DEFAULT and M_DEFAULT_VAL in both parameter searches.
#define M_DELTA_MEDIAN_THRESHOLD_DEFAULT 0.05_real64
#define M_DELTA_MAX_THRESHOLD_DEFAULT 0.10_real64
#define M_DELTA_EPSILON_DEFAULT 1.0e-10_real64
#define M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT 2_int32
