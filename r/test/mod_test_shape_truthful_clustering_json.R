source("r/load_tensor_omics.R")
source("r/test_helpers.R")

STOP_REASON_FIXED_POINT <- 4L

common_args <- list(
  k_min = 3L, k_density = 4L, chordal_dist_max_as_prcnt_of_range = 0.1, d_max = 1L,
  G_max = 2.0, RMSE_change_max = 0.5, f_max = 0.8, a = 3L,
  exclusion_radius_percentile = 50.0, bandwidth_percentile = 68.0,
  reconciliation_mode = "merge_overlap_coefficient", min_overlap_coefficient = 0.5
)

# D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which overlap
# on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble.
fixture <- function() {
  vectors <- matrix(c(0, 0, 1, 0, 2, 0, 3, 0), nrow = 2, ncol = 4)
  dim_names <- c("x", "y")
  seed_selection_mask <- c(TRUE, FALSE, FALSE, TRUE)

  ensemble_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_masks[1:3, 1] <- TRUE
  ensemble_masks[2:4, 2] <- TRUE

  ensemble_stop_reason <- c(STOP_REASON_FIXED_POINT, STOP_REASON_FIXED_POINT)
  ensemble_growth_radii <- c(1.0, 1.0)

  ensemble_k_history <- matrix(c(2, 3, 2, 3), nrow = 2, ncol = 2)
  ensemble_d_history <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  ensemble_G_history <- matrix(c(2.0, 1.5, 2.0, 1.5), nrow = 2, ncol = 2)

  ensemble_mu_history <- array(0.0, dim = c(2, 2, 2))
  ensemble_mu_history[, 1, 1] <- c(0.5, 0.0)
  ensemble_mu_history[, 2, 1] <- c(1.0, 0.0)
  ensemble_mu_history[, 1, 2] <- c(2.5, 0.0)
  ensemble_mu_history[, 2, 2] <- c(3.0, 0.0)

  ensemble_S_history <- array(0.0, dim = c(2, 2, 2))
  ensemble_S_history[1, 2, 1] <- 0.5

  ensemble_U_history <- array(0.0, dim = c(2, 2, 2, 2))
  ensemble_U_history[, 1, 2, 1] <- c(1.0, 0.0)

  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_low_confidence_masks[1, 1] <- TRUE

  super_ensembles <- matrix(c(1, 2, 0, 0), nrow = 2, ncol = 2)

  list(
    vectors = vectors, dim_names = dim_names, seed_selection_mask = seed_selection_mask,
    ensemble_masks = ensemble_masks, ensemble_stop_reason = ensemble_stop_reason,
    ensemble_growth_radii = ensemble_growth_radii, ensemble_U_history = ensemble_U_history,
    ensemble_S_history = ensemble_S_history, ensemble_d_history = ensemble_d_history,
    ensemble_G_history = ensemble_G_history, ensemble_mu_history = ensemble_mu_history,
    ensemble_k_history = ensemble_k_history, ensemble_low_confidence_masks = ensemble_low_confidence_masks,
    super_ensembles = super_ensembles
  )
}

call_serialize <- function(fx, filename, n_super_ensembles = 1L, ...) {
  args <- c(list(
    filename = filename, n_super_ensembles = n_super_ensembles, vectors = fx$vectors,
    dim_names = fx$dim_names, seed_selection_mask = fx$seed_selection_mask,
    ensemble_masks = fx$ensemble_masks, ensemble_stop_reason = fx$ensemble_stop_reason,
    ensemble_growth_radii = fx$ensemble_growth_radii, ensemble_U_history = fx$ensemble_U_history,
    ensemble_S_history = fx$ensemble_S_history, ensemble_d_history = fx$ensemble_d_history,
    ensemble_G_history = fx$ensemble_G_history, ensemble_mu_history = fx$ensemble_mu_history,
    ensemble_k_history = fx$ensemble_k_history, ensemble_low_confidence_masks = fx$ensemble_low_confidence_masks,
    super_ensembles = fx$super_ensembles
  ), common_args, list(...))
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

  assert_true(grepl(paste0(
    '"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],',
    '"n_ensembles":1,"n_low_confidence_ensembles":1,"ensembles":[1],',
    '"low_confidence_ensembles":[1],"seed_of":[1]'
  ), content, fixed = TRUE))
  assert_true(grepl(paste0(
    '"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],',
    '"n_ensembles":2,"n_low_confidence_ensembles":0,"ensembles":[1,2],',
    '"low_confidence_ensembles":[],"seed_of":[]'
  ), content, fixed = TRUE))

  assert_true(grepl(paste0(
    '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,',
    '"size":3,"d":1,"G":1.5000000000000000E+000,"mu":[1.0000000000000000E+000,0.0000000000000000E+000],',
    '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,',
    '"super_ensemble_id":1}'
  ), content, fixed = TRUE))
  assert_true(grepl(paste0(
    '{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,',
    '"size":3,"d":0,"G":1.5000000000000000E+000,"mu":[3.0000000000000000E+000,0.0000000000000000E+000],',
    '"super_ensemble_id":1}'
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
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  filename <- "test_stc_zero_ensembles_r.json"
  args <- c(list(
    filename = filename, n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks, super_ensembles = super_ensembles
  ), common_args)
  do.call(serialize_stc_results_as_json, args)

  content <- read_whole_file(filename)
  file.remove(filename)

  assert_true(grepl('"ensembles":[]', content, fixed = TRUE))
  assert_true(grepl('"super_ensembles":[]', content, fixed = TRUE))
  assert_true(grepl('"overlap_coefficient_matrix":[]', content, fixed = TRUE))
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
    ensemble_low_confidence_masks = fx$ensemble_low_confidence_masks, super_ensembles = fx$super_ensembles
  ), common_args)
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
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  args <- c(list(
    filename = "test_stc_invalid_r.json", n_super_ensembles = 0L, vectors = vectors, dim_names = dim_names,
    seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
    ensemble_stop_reason = ensemble_stop_reason, ensemble_growth_radii = ensemble_growth_radii,
    ensemble_U_history = ensemble_U_history, ensemble_S_history = ensemble_S_history,
    ensemble_d_history = ensemble_d_history, ensemble_G_history = ensemble_G_history,
    ensemble_mu_history = ensemble_mu_history, ensemble_k_history = ensemble_k_history,
    ensemble_low_confidence_masks = ensemble_low_confidence_masks, super_ensembles = super_ensembles
  ), common_args)

  assert_error(do.call(serialize_stc_results_as_json, args),
               "n_dimensions=0 must be rejected", ERR_EMPTY_INPUT)
}

run_all_tests()
