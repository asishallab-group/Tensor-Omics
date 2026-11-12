source("r/load_tensor_omics.R")

vectors <- matrix(c(0, 0.1, 0.2, 50, 55, 60), nrow = 1)
dimension_order <- 1L
kd_indices <- build_kd_index(vectors, dimension_order)
labels <- calculate_labels_as_density(vectors, 0.15, dimension_order, kd_indices, n_tiles = 2L)
stopifnot(all(labels == c(2, 3, 2, 1, 1, 1)))

# Hold the density labels and the seeds fixed to isolate the growth radius.
density_labels <- c(8, 8, 8, 4, 4, 4)
seed_indices <- c(1L, 4L)
global_result <- obtain_ensembles_global_radius(vectors, dimension_order, kd_indices, density_labels,
                                         seed_indices, r = 0.15, t_observables = 10L, n_tiles = 2L)
local_result <- obtain_ensembles_seed_radius(vectors, dimension_order, kd_indices, density_labels,
                                      seed_indices, seed_radii = c(0.15, 5.5), t_observables = 10L, n_tiles = 2L)
stopifnot(all(colSums(global_result$ensemble_matrix) == c(3, 1)),
          all(colSums(local_result$ensemble_matrix) == c(3, 3)),
          all(global_result$stop_reasons == c(1L, 4L)),
          all(local_result$stop_reasons == c(1L, 1L)),
          global_result$mad_ambient == 2, local_result$mad_ambient == 2,
          identical(merge_ensembles(global_result$ensemble_matrix),
                    global_result$ensemble_matrix[, 1, drop = FALSE]))

# The in/out logical mask returns a new R value and preserves the caller's input.
mask <- c(TRUE, rep(FALSE, 5))
grown <- grow_ensemble(vectors, mask, 0.15, dimension_order, kd_indices, density_labels)
stopifnot(identical(mask, c(TRUE, rep(FALSE, 5))),
          identical(grown, c(TRUE, TRUE, rep(FALSE, 4))))

history <- matrix(0, nrow = 4, ncol = 2)
bad_update <- tryCatch(compute_ensemble_observable(mask, density_labels, 0L, 1L, history), error = identity)
bad_accept <- tryCatch(accept_ensemble(history, 2L), error = identity)
stopifnot(inherits(bad_update, "error"), inherits(bad_accept, "error"),
          grepl("observables", conditionMessage(bad_update), fixed = TRUE),
          grepl("observables", conditionMessage(bad_accept), fixed = TRUE))

stopifnot(!("n_ensembles" %in% names(global_result)))
stopifnot(!("n_ensembles" %in% names(local_result)))
