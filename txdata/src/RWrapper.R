# tensoromics.R — R-Wrapper & Tests für libtensoromics.so

dyn.load("build/libtensoromics.so")


create_estimates <- function(n_tissues, n_genes, meta_n_rows, meta_max_char, meta_col_types) {
  list(
    n_tissues      = as.integer(n_tissues),
    n_genes        = as.integer(n_genes),
    meta_n_rows    = as.integer(meta_n_rows),
    meta_max_char  = as.integer(meta_max_char),
    meta_col_types = as.integer(meta_col_types),
    n_cols         = as.integer(length(meta_col_types)))
}

init_tensoromics <- function(est) {
  res <- .Fortran("init",
                  as.integer(est$n_tissues),
                  as.integer(est$n_genes),
                  as.integer(est$meta_n_rows),
                  as.integer(est$meta_max_char),
                  as.integer(est$meta_col_types),
                  as.integer(est$n_cols),
                  vec_out    = double(est$n_tissues * est$n_genes),
                  shift_out  = double(est$n_tissues * est$n_genes),
                  next_idx   = integer(1))
  
  list(
    vec_container = matrix(res$vec_out, nrow = est$n_tissues, ncol = est$n_genes),
    shift_vecs    = matrix(res$shift_out, nrow = est$n_tissues, ncol = est$n_genes),
    next_idx      = res$next_idx)
}

calculate_memory <- function(est) {
  .Fortran("calculate_memory_requirements",
           as.integer(est$n_tissues),
           as.integer(est$n_genes),
           as.integer(est$meta_n_rows),
           as.integer(est$meta_max_char),
           as.integer(est$meta_col_types),
           as.integer(est$n_cols),
           mem_bytes = integer(1))$mem_bytes
}

update_tensoromics <- function(tom, patch) {
  res <- .Fortran("update",
                  as.integer(nrow(tom$vec_container)),
                  as.integer(ncol(tom$vec_container)),
                  as.double(patch),
                  as.integer(ncol(patch)),
                  vec_container = as.double(tom$vec_container),
                  shift_vecs    = as.double(tom$shift_vecs),
                  next_idx      = as.integer(tom$next_idx),
                  indices       = integer(ncol(patch)))
  
  list(
    vec_container = matrix(res$vec_container, nrow = nrow(tom$vec_container)),
    shift_vecs    = matrix(res$shift_vecs, nrow = nrow(tom$vec_container)),
    next_idx      = res$next_idx,
    indices       = res$indices)
}

save_tensoromics <- function(tom, filename) {
  .Fortran("save",
           as.integer(nrow(tom$vec_container)),
           as.integer(ncol(tom$vec_container)),
           as.double(tom$vec_container),
           as.double(tom$shift_vecs),
           as.character(filename),
           as.integer(nchar(filename)))
}

# --------------------
# TEST FUNCTIONS:
# --------------------
test_init_mem <- function() {
  est <- create_estimates(8, 50, 20, 16, c(1,2,3))
  tom <- init_tensoromics(est)
  stopifnot(all(dim(tom$vec_container) == c(8,50)),
            tom$next_idx == 1L)
  cat("Init OK. Memory:", calculate_memory(est), "bytes\n")
}

test_update <- function() {
  est <- create_estimates(4, 30, 10, 12, rep(2,3))
  tom <- init_tensoromics(est)
  patch <- matrix(rnorm(4*5), nrow=4, ncol=5)
  tom2 <- update_tensoromics(tom, patch)
  stopifnot(tom2$next_idx == 6L)
  cat("Update OK. Indices:", tom2$indices, "\n")
}

test_save <- function() {
  est <- create_estimates(3, 10, 5, 8, rep(2,2))
  tom <- init_tensoromics(est)
  patch <- matrix(runif(3*2), nrow=3, ncol=2)
  tom2 <- update_tensoromics(tom, patch)
  fn <- file.path("data/test_tensoromics.txdata")
  save_tensoromics(tom2, fn)
  info <- file.info(fn)
  stopifnot(!is.na(info$size) && info$size > 0)
  cat("Save OK. File:", fn, "\n")
}

run_all_tests <- function() {
  cat("=== Testing init & memory ===\n"); test_init_mem()
  cat("=== Testing update       ===\n"); test_update()
  cat("=== Testing save         ===\n"); test_save()
  cat("=== All tests passed! ===\n")
}

if (interactive()) run_all_tests()

run_all_tests()
