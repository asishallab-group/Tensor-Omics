source("r/load_tensor_omics.R")
source("r/test_helpers.R")

STOP_REASON_FIXED_POINT <- 4L
STOP_REASON_REJECTED_AFTER_STABLE <- 2L

common_args <- list(
  k_min = 3L, k_density = 4L, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 1L,
  G_max = 2.0, RMSE_change_max = 0.5, f_max = 0.8, a = 3L,
  exclusion_radius_percentile = 50.0, bandwidth_percentile = 68.0,
  reconciliation_mode = "merge_overlap_coefficient", min_overlap_coefficient = 0.5
)

# D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which overlap
# on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble. Both ensembles
# are fully accepted at both retained history columns -- the rejected-trailing-column case gets
# its own dedicated fixture below. Every S_history entry that feeds an observable_history rmse
# is chosen so normal_error is a perfect square (4.0 or 9.0), avoiding any rounding-tie risk
# from the +epsilon guard around a normal_error of exactly 0.
fixture <- function() {
  vectors <- matrix(c(0, 0, 1, 0, 2, 0, 3, 0), nrow = 2, ncol = 4)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE, FALSE, FALSE, TRUE)

  ensemble_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_masks[1:3, 1] <- TRUE
  ensemble_masks[2:4, 2] <- TRUE

  ensemble_stop_reason <- c(STOP_REASON_FIXED_POINT, STOP_REASON_FIXED_POINT)
  ensemble_growth_radii <- c(1.0, 1.0)

  # k-1=1 at every retained column of both ensembles, to keep the eigenvalue arithmetic
  # (S**2/(k-1)) below simple; k_history's exact value is otherwise unasserted.
  ensemble_k_history <- matrix(c(2, 2, 2, 2), nrow = 2, ncol = 2)
  ensemble_d_history <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  ensemble_G_history <- matrix(c(2.0, 1.5, 2.0, 1.5), nrow = 2, ncol = 2)

  ensemble_mu_history <- array(0.0, dim = c(2, 2, 2))
  ensemble_mu_history[, 1, 1] <- c(0.5, 0.0)
  ensemble_mu_history[, 2, 1] <- c(1.0, 0.0)
  ensemble_mu_history[, 1, 2] <- c(2.5, 0.0)
  ensemble_mu_history[, 2, 2] <- c(3.0, 0.0)

  ensemble_S_history <- array(0.0, dim = c(2, 2, 2))
  ensemble_S_history[1, 1, 1] <- 2.0 # ensemble 1, iter 1 (d=0): normal_error=4.0, rmse=2.0
  ensemble_S_history[1, 2, 1] <- 0.5 # ensemble 1, iter 2 (d=1): s1, excluded from normal_error
  ensemble_S_history[2, 2, 1] <- 3.0 # ensemble 1, iter 2 (d=1): normal_error=9.0, rmse=3.0
  ensemble_S_history[1, 1, 2] <- 2.0 # ensemble 2, iter 1 (d=0): normal_error=4.0, rmse=2.0
  ensemble_S_history[1, 2, 2] <- 3.0 # ensemble 2, iter 2 (d=0): normal_error=9.0, rmse=3.0

  ensemble_U_history <- array(0.0, dim = c(2, 2, 2, 2))
  ensemble_U_history[, 1, 2, 1] <- c(1.0, 0.0) # ensemble 1's u1 (tangent)
  ensemble_U_history[, 2, 2, 1] <- c(0.0, 1.0) # ensemble 1's normal direction (residual_length)
  ensemble_U_history[, 1, 2, 2] <- c(1.0, 0.0) # ensemble 2 has d=0: only ever read by residual_length
  ensemble_U_history[, 2, 2, 2] <- c(0.0, 1.0)

  # Both ensembles accepted at every retained column -- no rejection in this fixture.
  ensemble_accepted_history <- matrix(TRUE, nrow = 2, ncol = 2)

  # Ensemble 1: seed = point 1; point 2 joins at iteration 1, point 3 at iteration 2; point 4
  # never joins. T (max) = 2. Ensemble 2: seed = point 4; point 3 joins at iteration 1, point 2
  # at iteration 2; point 1 never joins. T (max) = 2.
  ensemble_member_added_at_step <- matrix(c(0, 1, 2, -1, -1, 2, 1, 0), nrow = 4, ncol = 2)

  # Bootstrap (iteration 1) basis, duplicating history column 1 exactly -- both d=0 here, so
  # U_first is never actually read (stc_chordal_distance requires d>0 to be "applicable").
  ensemble_U_first <- array(0.0, dim = c(2, 2, 2))
  ensemble_U_first[, , 1] <- ensemble_U_history[, , 1, 1]
  ensemble_U_first[, , 2] <- ensemble_U_history[, , 1, 2]
  ensemble_d_first <- c(ensemble_d_history[1, 1], ensemble_d_history[1, 2])

  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_low_confidence_masks[1, 1] <- TRUE

  super_ensembles <- matrix(c(1, 2, 0, 0), nrow = 2, ncol = 2)

  list(
    vectors = vectors, dim_names = dim_names, seed_selection_mask = seed_selection_mask,
    ensemble_masks = ensemble_masks, ensemble_stop_reason = ensemble_stop_reason,
    ensemble_growth_radii = ensemble_growth_radii, ensemble_U_history = ensemble_U_history,
    ensemble_S_history = ensemble_S_history, ensemble_d_history = ensemble_d_history,
    ensemble_G_history = ensemble_G_history, ensemble_mu_history = ensemble_mu_history,
    ensemble_k_history = ensemble_k_history, ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  )
}

# The new required per-ensemble reconciliation-eligibility arrays, all-eligible unless a test
# overrides them via ... -- mirrors test/mod_test_shape_truthful_clustering_json.F90's own
# per-test `ensemble_eligible = .true.` convention.
all_eligible_args <- function(n) {
  list(
    ensemble_eligible = rep(TRUE, n), ensemble_eligible_by_stop_condition = rep(TRUE, n),
    ensemble_eligible_by_dimension = rep(TRUE, n), ensemble_eligible_by_var_explained = rep(TRUE, n)
  )
}

call_serialize <- function(fx, filename, n_super_ensembles = 1L, ...) {
  overrides <- list(...)
  n_selected_seed <- ncol(fx$ensemble_masks)
  eligible <- all_eligible_args(n_selected_seed)
  eligible[names(overrides)[names(overrides) %in% names(eligible)]] <- NULL
  args <- c(list(
    filename = filename, n_super_ensembles = n_super_ensembles, vectors = fx$vectors,
    dim_names = fx$dim_names, seed_selection_mask = fx$seed_selection_mask,
    ensemble_masks = fx$ensemble_masks, ensemble_stop_reason = fx$ensemble_stop_reason,
    ensemble_growth_radii = fx$ensemble_growth_radii, ensemble_U_history = fx$ensemble_U_history,
    ensemble_S_history = fx$ensemble_S_history, ensemble_d_history = fx$ensemble_d_history,
    ensemble_G_history = fx$ensemble_G_history, ensemble_mu_history = fx$ensemble_mu_history,
    ensemble_k_history = fx$ensemble_k_history, ensemble_accepted_history = fx$ensemble_accepted_history,
    ensemble_member_added_at_step = fx$ensemble_member_added_at_step,
    ensemble_low_confidence_masks = fx$ensemble_low_confidence_masks,
    ensemble_U_first = fx$ensemble_U_first, ensemble_d_first = fx$ensemble_d_first,
    super_ensembles = fx$super_ensembles
  ), common_args, eligible, overrides)
  do.call(serialize_stc_results_as_json, args)
}

read_whole_file <- function(filename) {
  paste(readLines(filename, warn = FALSE), collapse = "\n")
}

test_json_two_ensembles_with_overlap <- function() {
  fx <- fixture()
  filename <- "test_stc_two_ensembles_r.json"
  call_serialize(fx, filename)
  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"dim_names":["x","y"]', content, fixed = TRUE))
  assert_true(grepl('"n_vectors":4', content, fixed = TRUE))
  assert_true(grepl('"n_dimensions":2', content, fixed = TRUE))
  assert_true(grepl('"n_ensembles":2', content, fixed = TRUE))
  assert_true(grepl('"k_min":3', content, fixed = TRUE))
  assert_true(grepl('"reconciliation_mode":"merge_overlap_coefficient"', content, fixed = TRUE))

  # residual_length is the point's distance off each ensemble's *normal* subspace at its final
  # retained basis (identity here, see fixture()); point 1 sits exactly at ensemble 1's mean
  # (residual 0); point 2 sits at ensemble 1's mean too (residual 0) but 2.0 off ensemble 2's
  # mean along its own (identity) normal directions.
  assert_true(grepl(paste0(
    '"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],',
    '"n_ensembles":1,"n_low_confidence_ensembles":1,',
    '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}],',
    '"low_confidence_ensembles":[1],"seed_of":[1]'
  ), content, fixed = TRUE))
  assert_true(grepl(paste0(
    '"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],',
    '"n_ensembles":2,"n_low_confidence_ensembles":0,',
    '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000},',
    '{"id":2,"residual_length":2.0000000000000000E+000}],',
    '"low_confidence_ensembles":[],"seed_of":[]'
  ), content, fixed = TRUE))

  # ensembles: #1 has d=1 (u1/s1/line_start/line_end present), #2 has d=0 (no tangent direction
  # or line at all). Neither ensemble's drift is ever computable here (ensemble 1's
  # iter1-vs-iter2 d mismatches, 0 vs 1; ensemble 2's two iterations are both d=0, and
  # stc_chordal_distance requires d>0 to be "applicable" at all).
  assert_true(grepl(paste0(
    '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,',
    '"size":3,"t_final":2,"d":1,"G":1.5000000000000000E+000,',
    '"mu":[1.0000000000000000E+000,0.0000000000000000E+000],',
    '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,',
    '"line_start":[0.0000000000000000E+000,0.0000000000000000E+000],',
    '"line_end":[2.0000000000000000E+000,0.0000000000000000E+000],',
    '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,',
    '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,',
    '"rmse":3.0000000000000000E+000,"drift":null}],',
    '"super_ensemble_id":1,"final_chordal_distance":null,',
    '"reconciliation_eligible":true,"excluded_by":[]}'
  ), content, fixed = TRUE))
  assert_true(grepl(paste0(
    '{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,',
    '"size":3,"t_final":2,"d":0,"G":1.5000000000000000E+000,',
    '"mu":[3.0000000000000000E+000,0.0000000000000000E+000],',
    '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,',
    '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,',
    '"rmse":3.0000000000000000E+000,"drift":null}],',
    '"super_ensemble_id":1,"final_chordal_distance":null,',
    '"reconciliation_eligible":true,"excluded_by":[]}'
  ), content, fixed = TRUE))
  assert_true(!grepl('"u2"', content, fixed = TRUE))

  assert_true(grepl('"super_ensembles":[{"group_id":1,"ensemble_ids":[1,2]}]', content, fixed = TRUE))
  assert_true(grepl(
    '"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":6.6666666666666663E-001}]',
    content, fixed = TRUE
  ))
}

test_json_estimated_params_included <- function() {
  fx <- fixture()
  filename <- "test_stc_estimated_params_r.json"
  call_serialize(fx, filename,
                 estimated_k_min = 5L, estimated_k_density = 6L, estimated_density_quantile = 0.75,
                 estimated_chordal_dist_max_as_prcnt_of_range = 0.2, estimated_G_max = 3.0, estimated_d_max = 2L)
  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"estimated_k_min":5', content, fixed = TRUE))
  assert_true(grepl('"estimated_k_density":6', content, fixed = TRUE))
  assert_true(grepl('"estimated_d_max":2', content, fixed = TRUE))
}

# n_selected_seed=0, derived from an all-FALSE seed_selection_mask, is a valid, well-defined
# "no ensembles" input, not an error.
test_json_zero_ensembles <- function() {
  vectors <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(FALSE, FALSE)
  ensemble_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  ensemble_stop_reason <- integer(0)
  ensemble_growth_radii <- numeric(0)
  ensemble_U_history <- array(0.0, dim = c(2, 2, 1, 0))
  ensemble_S_history <- array(0.0, dim = c(2, 1, 0))
  ensemble_d_history <- array(0L, dim = c(1, 0))
  ensemble_G_history <- array(0.0, dim = c(1, 0))
  ensemble_mu_history <- array(0.0, dim = c(2, 1, 0))
  ensemble_k_history <- array(0L, dim = c(1, 0))
  ensemble_accepted_history <- array(FALSE, dim = c(1, 0))
  ensemble_member_added_at_step <- array(0L, dim = c(2, 0))
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  ensemble_U_first <- array(0.0, dim = c(2, 2, 0))
  ensemble_d_first <- integer(0)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  filename <- "test_stc_zero_ensembles_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(0))
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"ensembles":[]', content, fixed = TRUE))
  assert_true(grepl('"super_ensembles":[]', content, fixed = TRUE))
  assert_true(grepl('"overlap_coefficient_matrix":[]', content, fixed = TRUE))
}

# An ensemble whose k_history is entirely zero never produced an observable at all (only
# possible for STOP_REASON_MAX_SIZE firing at the bootstrap step itself). Its
# t_final/d/G/mu/tangent/observable_history keys must be omitted, not merely null --
# final_chordal_distance is the one exception, always present (null here). The one point that
# is a member of this ensemble still gets its usual `ensembles` entry, with a residual_length
# of 0.0 (the documented fallback when no retained basis exists at all).
test_json_no_history_ensemble_omits_observable_keys <- function() {
  vectors <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE, FALSE)
  ensemble_masks <- matrix(c(TRUE, FALSE), nrow = 2, ncol = 1)
  ensemble_stop_reason <- c(STOP_REASON_FIXED_POINT)
  ensemble_growth_radii <- c(1.0)
  ensemble_k_history <- matrix(0L, nrow = 1, ncol = 1)
  ensemble_d_history <- array(0L, dim = c(1, 1))
  ensemble_G_history <- array(0.0, dim = c(1, 1))
  ensemble_mu_history <- array(0.0, dim = c(2, 1, 1))
  ensemble_S_history <- array(0.0, dim = c(2, 1, 1))
  ensemble_U_history <- array(0.0, dim = c(2, 2, 1, 1))
  ensemble_accepted_history <- array(TRUE, dim = c(1, 1))
  ensemble_member_added_at_step <- matrix(c(0L, -1L), nrow = 2, ncol = 1)
  ensemble_U_first <- array(0.0, dim = c(2, 2, 1))
  ensemble_d_first <- c(0L)
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 1)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  filename <- "test_stc_no_history_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(1))
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl(paste0(
    '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,',
    '"size":1,"super_ensemble_id":null,"final_chordal_distance":null,',
    '"reconciliation_eligible":true,"excluded_by":[]}'
  ), content, fixed = TRUE))
  assert_true(!grepl('"t_final"', content, fixed = TRUE))
  assert_true(!grepl('"d":', content, fixed = TRUE))
  assert_true(!grepl('"G":', content, fixed = TRUE))
  assert_true(!grepl('"mu":', content, fixed = TRUE))
  assert_true(!grepl('"line_start"', content, fixed = TRUE))
  assert_true(!grepl('"observable_history"', content, fixed = TRUE))
  assert_true(grepl('"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}]', content, fixed = TRUE))
}

test_html_report_wraps_json_in_template_and_d3 <- function() {
  fx <- fixture()
  filename <- "test_stc_report_r.html"
  args <- c(list(
    filename = filename, n_super_ensembles = 1L, vectors = fx$vectors, dim_names = fx$dim_names,
    seed_selection_mask = fx$seed_selection_mask, ensemble_masks = fx$ensemble_masks,
    ensemble_stop_reason = fx$ensemble_stop_reason, ensemble_growth_radii = fx$ensemble_growth_radii,
    ensemble_U_history = fx$ensemble_U_history, ensemble_S_history = fx$ensemble_S_history,
    ensemble_d_history = fx$ensemble_d_history, ensemble_G_history = fx$ensemble_G_history,
    ensemble_mu_history = fx$ensemble_mu_history, ensemble_k_history = fx$ensemble_k_history,
    ensemble_accepted_history = fx$ensemble_accepted_history,
    ensemble_member_added_at_step = fx$ensemble_member_added_at_step,
    ensemble_low_confidence_masks = fx$ensemble_low_confidence_masks,
    ensemble_U_first = fx$ensemble_U_first, ensemble_d_first = fx$ensemble_d_first,
    super_ensembles = fx$super_ensembles
  ), common_args, all_eligible_args(2))
  do.call(write_stc_interactive_html_report, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(startsWith(content, "<!DOCTYPE html>"))
  assert_true(grepl("d3js.org", content, fixed = TRUE))
  assert_true(grepl('const DATA = {"dim_names":["x","y"]', content, fixed = TRUE))
  assert_true(grepl("</html>", content, fixed = TRUE))
}

test_json_invalid_n_dimensions <- function() {
  vectors <- matrix(numeric(0), nrow = 0, ncol = 2)
  dim_names <- character(0)
  seed_selection_mask <- c(FALSE, FALSE)
  ensemble_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  ensemble_stop_reason <- integer(0)
  ensemble_growth_radii <- numeric(0)
  ensemble_U_history <- array(0.0, dim = c(0, 0, 1, 0))
  ensemble_S_history <- array(0.0, dim = c(0, 1, 0))
  ensemble_d_history <- array(0L, dim = c(1, 0))
  ensemble_G_history <- array(0.0, dim = c(1, 0))
  ensemble_mu_history <- array(0.0, dim = c(0, 1, 0))
  ensemble_k_history <- array(0L, dim = c(1, 0))
  ensemble_accepted_history <- array(FALSE, dim = c(1, 0))
  ensemble_member_added_at_step <- array(0L, dim = c(2, 0))
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  ensemble_U_first <- array(0.0, dim = c(0, 0, 0))
  ensemble_d_first <- integer(0)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  args <- c(list(
    filename = "test_stc_invalid_r.json", n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(0))

  assert_error(do.call(serialize_stc_results_as_json, args),
               "n_dimensions=0 must be rejected", ERR_EMPTY_INPUT)
}

# Regression test for a real, pre-existing bug: stc_push_ensemble_history also pushes a
# *rejected* final candidate into the trailing history window right before growth halts on
# STOP_REASON_REJECTED_IMMEDIATELY/STOP_REASON_REJECTED_AFTER_STABLE -- so the last *populated*
# history column is not always the ensemble's real last *accepted* state. This fixture's
# column 1 is the true accepted state (d=1, G=1.0, mu=[0.5,0.0], u1=[1,0]); its column 2 is a
# rejected candidate with deliberately different, wrong-if-leaked values (d=0, G=99.0,
# mu=[9.0,9.0], accepted_history=FALSE). Every reported "current state" field must reflect
# column 1, and observable_history must contain exactly that one entry.
test_json_rejected_trailing_column_uses_last_accepted <- function() {
  vectors <- matrix(c(0, 0, 1, 0, 9, 9), nrow = 2, ncol = 3)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE, FALSE, FALSE)
  ensemble_masks <- matrix(c(TRUE, TRUE, FALSE), nrow = 3, ncol = 1)
  ensemble_stop_reason <- c(STOP_REASON_REJECTED_AFTER_STABLE)
  ensemble_growth_radii <- c(1.0)

  ensemble_k_history <- matrix(c(2, 3), nrow = 2, ncol = 1)
  ensemble_d_history <- matrix(c(1, 0), nrow = 2, ncol = 1)
  ensemble_G_history <- matrix(c(1.0, 99.0), nrow = 2, ncol = 1)
  ensemble_mu_history <- array(0.0, dim = c(2, 2, 1))
  ensemble_mu_history[, 1, 1] <- c(0.5, 0.0)
  ensemble_mu_history[, 2, 1] <- c(9.0, 9.0)
  ensemble_S_history <- array(0.0, dim = c(2, 2, 1))
  ensemble_S_history[, 1, 1] <- c(0.5, 2.0) # eigen(2) = 4.0/(2-1) = 4.0 -> rmse = 2.0
  ensemble_S_history[, 2, 1] <- c(7.0, 7.0)
  ensemble_U_history <- array(0.0, dim = c(2, 2, 2, 1))
  ensemble_U_history[, 1, 1, 1] <- c(1.0, 0.0)
  ensemble_U_history[, 2, 1, 1] <- c(0.0, 1.0)
  ensemble_accepted_history <- matrix(c(TRUE, FALSE), nrow = 2, ncol = 1)

  # T=1: only the bootstrap was ever actually accepted.
  ensemble_member_added_at_step <- matrix(c(0, 1, -1), nrow = 3, ncol = 1)

  ensemble_U_first <- array(0.0, dim = c(2, 2, 1))
  ensemble_U_first[, , 1] <- ensemble_U_history[, , 1, 1]
  ensemble_d_first <- c(1L)
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 3, ncol = 1)
  super_ensembles <- matrix(0L, nrow = 1, ncol = 0)

  filename <- "test_stc_rejected_trailing_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(1))
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"stop_reason":"rejected_after_stable"', content, fixed = TRUE))
  assert_true(grepl('"size":2', content, fixed = TRUE))
  assert_true(grepl('"t_final":1', content, fixed = TRUE)) # the rejected step never counted
  assert_true(grepl('"d":1', content, fixed = TRUE)) # from the accepted column, not the rejected column's 0
  assert_true(grepl('"G":1.0000000000000000E+000', content, fixed = TRUE)) # not the rejected column's 99.0
  assert_true(grepl('"mu":[5.0000000000000000E-001,0.0000000000000000E+000]', content, fixed = TRUE)) # not [9,9]
  assert_true(grepl('"u1":[1.0000000000000000E+000,0.0000000000000000E+000]', content, fixed = TRUE))
  assert_true(grepl(paste0(
    '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,',
    '"rmse":2.0000000000000000E+000,"drift":null}]'
  ), content, fixed = TRUE))
  assert_true(!grepl("9.9000000000000000E+001", content, fixed = TRUE)) # rejected G=99.0 never appears anywhere
  # Point 3's own coords legitimately contain [9,9], so check the ensemble's mu specifically.
  assert_true(!grepl('"mu":[9.0000000000000000E+000,9.0000000000000000E+000]', content, fixed = TRUE))
}

# A 3-retained-iteration, d=1-throughout ensemble whose tangent direction rotates 90 degrees,
# then back: u1 = [1,0] (iter 1 = bootstrap = U_first), [0,1] (iter 2), [1,0] (iter 3). Every
# consecutive pair is orthogonal, so the chordal distance sqrt(1 - cos(theta)**2) is exactly
# 1.0 for both drift entries. final_chordal_distance compares iteration 3 against
# {U_first, iter1, iter2}: 0.0, 0.0, 1.0 -- max = 1.0.
test_json_drift_and_final_chordal_distance <- function() {
  vectors <- matrix(c(0, 0, 5, 5), nrow = 2, ncol = 2)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE, FALSE)
  ensemble_masks <- matrix(c(TRUE, TRUE), nrow = 2, ncol = 1)
  ensemble_stop_reason <- c(STOP_REASON_FIXED_POINT)
  ensemble_growth_radii <- c(1.0)

  ensemble_k_history <- matrix(c(2, 2, 2), nrow = 3, ncol = 1)
  ensemble_d_history <- matrix(c(1, 1, 1), nrow = 3, ncol = 1)
  ensemble_G_history <- matrix(c(1.0, 2.0, 3.0), nrow = 3, ncol = 1)
  ensemble_mu_history <- array(0.0, dim = c(2, 3, 1))

  # s1=1.0 throughout (unasserted); the second singular value drives rmse via eigen(2) (d=1
  # excludes eigen(1)): 2**2/1=4 -> rmse=2, 3**2/1=9 -> rmse=3, 4**2/1=16 -> rmse=4.
  ensemble_S_history <- array(0.0, dim = c(2, 3, 1))
  ensemble_S_history[, 1, 1] <- c(1.0, 2.0)
  ensemble_S_history[, 2, 1] <- c(1.0, 3.0)
  ensemble_S_history[, 3, 1] <- c(1.0, 4.0)

  ensemble_U_history <- array(0.0, dim = c(2, 2, 3, 1))
  ensemble_U_history[, 1, 1, 1] <- c(1.0, 0.0)
  ensemble_U_history[, 1, 2, 1] <- c(0.0, 1.0)
  ensemble_U_history[, 1, 3, 1] <- c(1.0, 0.0)

  ensemble_accepted_history <- array(TRUE, dim = c(3, 1))
  ensemble_member_added_at_step <- matrix(c(0, 3), nrow = 2, ncol = 1) # T=3

  ensemble_U_first <- array(0.0, dim = c(2, 2, 1))
  ensemble_U_first[, 1, 1] <- c(1.0, 0.0) # matches iteration 1 exactly
  ensemble_d_first <- c(1L)
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 1)
  super_ensembles <- matrix(0L, nrow = 1, ncol = 0)

  filename <- "test_stc_drift_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(1))
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"t_final":3', content, fixed = TRUE))
  assert_true(grepl(paste0(
    '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,',
    '"rmse":2.0000000000000000E+000,"drift":null},',
    '{"iteration":2,"g":2.0000000000000000E+000,"rmse":3.0000000000000000E+000,',
    '"drift":1.0000000000000000E+000},',
    '{"iteration":3,"g":3.0000000000000000E+000,"rmse":4.0000000000000000E+000,',
    '"drift":1.0000000000000000E+000}]'
  ), content, fixed = TRUE))
  assert_true(grepl('"final_chordal_distance":1.0000000000000000E+000', content, fixed = TRUE))
}

# observable_history's rmse key is present but null when that column's ensemble size was <=1
# (k_history entry of 1 -- normal_error/(k-1) would divide by zero).
test_json_observable_history_rmse_null_when_size_le_one <- function() {
  vectors <- matrix(c(0, 0), nrow = 2, ncol = 1)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE)
  ensemble_masks <- matrix(c(TRUE), nrow = 1, ncol = 1)
  ensemble_stop_reason <- c(STOP_REASON_FIXED_POINT)
  ensemble_growth_radii <- c(1.0)
  ensemble_k_history <- matrix(1L, nrow = 1, ncol = 1) # degenerate: k-1=0
  ensemble_d_history <- array(0L, dim = c(1, 1))
  ensemble_G_history <- matrix(1.0, nrow = 1, ncol = 1)
  ensemble_mu_history <- array(0.0, dim = c(2, 1, 1))
  ensemble_S_history <- array(1.0, dim = c(2, 1, 1))
  ensemble_U_history <- array(0.0, dim = c(2, 2, 1, 1))
  ensemble_accepted_history <- array(TRUE, dim = c(1, 1))
  # T=1, so idx=1 maps to iteration 1, not the seed sentinel 0.
  ensemble_member_added_at_step <- matrix(1L, nrow = 1, ncol = 1)
  ensemble_U_first <- array(0.0, dim = c(2, 2, 1))
  ensemble_d_first <- c(0L)
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 1, ncol = 1)
  super_ensembles <- matrix(0L, nrow = 1, ncol = 0)

  filename <- "test_stc_rmse_null_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_accepted_history = ensemble_accepted_history,
    ensemble_member_added_at_step = ensemble_member_added_at_step,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    ensemble_U_first = ensemble_U_first, ensemble_d_first = ensemble_d_first,
    super_ensembles = super_ensembles
  ), common_args, all_eligible_args(1))
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl(paste0(
    '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,',
    '"rmse":null,"drift":null}]'
  ), content, fixed = TRUE))
}

# Reuses fixture()'s two intersecting ensembles (OC = 2/3), but overrides ensemble 2's stop
# reason to STOP_REASON_REJECTED_AFTER_STABLE. Without allowed_stop_reasons, both ensembles are
# eligible, the (1,2) pair appears in overlap_coefficient_matrix as usual, and
# params.excluded_stop_reasons is empty. With allowed_stop_reasons excluding
# STOP_REASON_REJECTED_AFTER_STABLE, ensemble 2 is eligible=FALSE (this module honors the
# ensemble_eligible* arguments directly, not derived internally from allowed_stop_reasons --
# that argument is reported for transparency only, in params.excluded_stop_reasons), so the
# (1,2) pair must vanish from overlap_coefficient_matrix entirely.
test_json_stop_reason_filter_excludes_pair <- function() {
  fx <- fixture()
  fx$ensemble_stop_reason[2] <- STOP_REASON_REJECTED_AFTER_STABLE

  filename <- "test_stc_filter_baseline_r.json"
  call_serialize(fx, filename)
  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"excluded_stop_reasons":[]', content, fixed = TRUE))
  assert_true(grepl(paste0(
    '"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":',
    '6.6666666666666663E-001}]'
  ), content, fixed = TRUE))
  assert_true(grepl('"reconciliation_eligible":true,"excluded_by":[]', content, fixed = TRUE))

  allowed <- c(TRUE, FALSE, TRUE, TRUE)
  filename <- "test_stc_filter_excluded_r.json"
  call_serialize(fx, filename, allowed_stop_reasons = allowed,
                 ensemble_eligible = c(TRUE, FALSE), ensemble_eligible_by_stop_condition = c(TRUE, FALSE))
  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"excluded_stop_reasons":["rejected_after_stable"]', content, fixed = TRUE))
  assert_true(grepl('"overlap_coefficient_matrix":[]', content, fixed = TRUE))
  assert_true(grepl('"id":2,"seed_point_id":4,"stop_reason":"rejected_after_stable"', content, fixed = TRUE))
}

# Exact-string-match regression test: ensemble 1 is fully eligible (excluded_by an empty
# array); ensemble 2 is excluded by two criteria at once (stop condition and dimension, not
# variance explained). Also proves the JSON key is genuinely `excluded_by` (11 chars), never
# the 27-char `reconciliation_excluded_by` that silently truncated in ensemble_keys's
# character(len=24) buffer during this module's own development.
test_json_reconciliation_eligible_and_excluded_by_exact_strings <- function() {
  fx <- fixture()

  filename <- "test_stc_eligible_excluded_by_r.json"
  call_serialize(fx, filename,
                 ensemble_eligible = c(TRUE, FALSE), ensemble_eligible_by_stop_condition = c(TRUE, FALSE),
                 ensemble_eligible_by_dimension = c(TRUE, FALSE), ensemble_eligible_by_var_explained = c(TRUE, TRUE))
  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl(paste0(
    '"super_ensemble_id":1,"final_chordal_distance":null,',
    '"reconciliation_eligible":true,"excluded_by":[]}'
  ), content, fixed = TRUE))
  assert_true(grepl('"reconciliation_eligible":false,"excluded_by":["stop_condition","dimension"]',
                    content, fixed = TRUE))
  assert_true(!grepl("reconciliation_excluded_by", content, fixed = TRUE))
}

run_all_tests()
