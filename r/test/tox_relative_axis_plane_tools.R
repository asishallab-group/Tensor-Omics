# Load shared library
source("r/tensoromics_functions.R")


test_omics_vector_RAP_projection_call <- function() {
  n_axes <- 5
  n_selected_axes <- 2
  n_vecs <- 10
  n_selected_vecs <- 5

  vecs <- matrix(runif(n_axes * n_vecs), nrow = n_axes)

  vecs_selection_mask <- integer(n_vecs)
  vecs_selection_mask[sample(n_vecs, n_selected_vecs)] <- 1
  n_selected_vecs <- sum(vecs_selection_mask)

  axes_selection_mask <- integer(n_axes)
  axes_selection_mask[sample(n_axes, n_selected_axes)] <- 1

  projections <- omics_vector_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

  for (i_vec in seq_len(n_selected_vecs)) {
    col_sum <- sum(projections[, i_vec])
    stopifnot(abs(col_sum) < 1e-6)
  }
}

test_omics_field_RAP_projection_call <- function() {
  n_axes <- 5
  n_selected_axes <- 2
  n_vecs <- 10
  n_selected_vecs <- 5

  vecs <- matrix(runif(2 * n_axes * n_vecs), nrow = 2 * n_axes)

  vecs_selection_mask <- integer(n_vecs)
  vecs_selection_mask[sample(n_vecs, n_selected_vecs)] <- 1
  n_selected_vecs <- sum(vecs_selection_mask)

  axes_selection_mask <- integer(n_axes)
  axes_selection_mask[sample(n_axes, n_selected_axes)] <- 1

  projections <- omics_field_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

  for (i_vec in seq_len(n_selected_vecs)) {
    col_sum <- sum(projections[, i_vec])
    stopifnot(abs(col_sum) < 1e-6)
  }
}

run_relative_axis_plane_tests <- function() {
  passed <- 0
  failed <- 0
  test_names <- c("omics_vector_RAP_projection_call", "omics_field_RAP_projection_call")
  test_functions <- list(test_omics_vector_RAP_projection_call, test_omics_field_RAP_projection_call)
  cat("Running R relative axis plane projection tests...\n")
  cat("============================================\n")
  for (i in seq_along(test_functions)) {
    cat(sprintf("\n--- %s ---\n", test_names[i]))
    tryCatch({
      test_functions[[i]]()
      cat(sprintf("✓ %s\n", test_names[i]))
      passed <- passed + 1
    }, error = function(e) {
      cat(sprintf("✗ FAILED: %s\n", e$message))
      failed <- failed + 1
    })
  }
  cat("\n============================================\n")
  cat(sprintf("Tests completed: %d passed, %d failed\n", passed, failed))
  if (failed == 0) {
    cat("🎉 All R relative axis plane projection tests passed!\n")
  } else {
    cat("❌ Some tests failed. Please check the output above.\n")
  }
}

run_relative_axis_plane_tests()
