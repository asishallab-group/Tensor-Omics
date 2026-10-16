! CM_ values that more than one module of the tox_data_integration family needs. A dummy extent or
! a DM_DEFAULT cannot name a Fortran parameter -- the generated C and R bindings need the number
! itself -- so each such value is a macro, written once, here; the module that owns it defines its
! parameter from the macro, and every module that needs the literal includes this file. A CM_ value
! only one module uses stays in that module.

! Fixed ceiling on any histogram bin count: MAX_N_BINS of tox_data_integration_js_comp_test_binning_impl,
! and the extent of the bin-sized work arrays of the bin-count, search and final-run modules.
#define CM_MAX_N_BINS 256

! `bin_count_mode`'s values (the MODE_BIN_COUNT_* of tox_data_integration_js_comp_test_binning_impl):
! how a reference point's histogram bin count is chosen. Occupancy is the default of all four
! pipeline entry points; determine_bin_count_stone_impl rejects it.
#define CM_MODE_BIN_COUNT_OCCUPANCY 0_int32
#define CM_MODE_BIN_COUNT_STONE_EQUAL 1_int32
#define CM_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED 2_int32
#define CM_MODE_BIN_COUNT_STONE_POOLED 3_int32

! Defaults of the occupancy bin-count search (determine_bin_count_occupancy_impl, also forwarded by
! the search and final-run modules). CM_OCCUPANCY_GAMMA_DEFAULT (growth of a candidate bin count)
! is unrelated to GAMMA (decay of the fixed-k grid's n_points).
#define CM_OCCUPANCY_M_MIN_DEFAULT 3_int32
#define CM_OCCUPANCY_M_MAX_DEFAULT 120_int32
#define CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT 10_int32
#define CM_OCCUPANCY_GAMMA_DEFAULT 1.25_real64
#define CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.05_real64
#define CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.95_real64

! Defaults of construct_adaptive_neighborhoods_impl's optional knobs (tox_data_integration_preprocessing_impl),
! shared by DM_DEFAULT and M_DEFAULT_VAL there and in the adaptive parameter search.
#define CM_ADAPTIVE_TAU_DEFAULT 0.1_real64
#define CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT 1.0_real64
#define CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT 0_int32
#define CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT 1_int32
