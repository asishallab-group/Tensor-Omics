# tensoromics.R — R-Wrapper & Tests für libtensoromics.so

dyn.load("libtensoromics.so")

create_estimates <- function(n_tissues, n_genes, meta_n_rows, meta_max_char, meta_col_types) {
  list(
    n_tissues      = as.integer(n_tissues),
    n_genes        = as.integer(n_genes),
    meta_n_rows    = as.integer(meta_n_rows),
    meta_max_char  = as.integer(meta_max_char),
    meta_col_types = as.integer(meta_col_types),
    n_cols         = as.integer(length(meta_col_types))
  )
}

init_tensoromics <- function(est) {
  res <- .Fortran("init_",
                  est$n_tissues, est$n_genes,
                  est$meta_n_rows, est$meta_max_char,
                  est$meta_col_types, est$n_cols,
                  vec_out    = double(est$n_tissues * est$n_genes),
                  shift_out  = double(est$n_tissues * est$n_genes),
                  next_idx   = integer(1))
  list(
    vec_container = matrix(res$vec_out,   nrow = est$n_tissues, ncol = est$n_genes),
    shift_vecs    = matrix(res$shift_out, nrow = est$n_tissues, ncol = est$n_genes),
    next_idx      = res$next_idx
  )
}

calculate_memory <- function(est) {
  .Fortran("calculate_memory_requirements_",
           est$n_tissues, est$n_genes,
           est$meta_n_rows, est$meta_max_char,
           est$meta_col_types, est$n_cols,
           mem_bytes = integer(1))$mem_bytes
}

update_tensoromics <- function(tom, patch) {
  n_tissues <- nrow(patch)
  n_genes   <- ncol(tom$vec_container)
  n_patch   <- ncol(patch)
  res <- .Fortran("update_",
                  as.integer(n_tissues),
                  as.integer(n_genes),
                  as.double(as.vector(patch)),
                  as.integer(n_patch),
                  as.double(as.vector(tom$vec_container)),
                  as.double(as.vector(tom$shift_vecs)),
                  as.integer(tom$next_idx),
                  indices = integer(n_patch))
  list(
    vec_container = matrix(res[[5]], nrow = n_tissues, ncol = n_genes),
    shift_vecs    = matrix(res[[6]], nrow = n_tissues, ncol = n_genes),
    next_idx      = res[[7]],
    indices       = res$indices
  )
}

save_tensoromics <- function(tom, filename) {
  .Fortran("save_",
           as.integer(nrow(tom$vec_container)),
           as.integer(ncol(tom$vec_container)),
           as.double(as.vector(tom$vec_container)),
           as.double(as.vector(tom$shift_vecs)),
           as.character(filename),
           as.integer(nchar(filename)))
}

# --------------------
# TEST-FUNKTIONEN:
# -------------------- 
test_init_mem <- function() {
  est <- create_estimates(8, 50, 20, 16, c(1,2,3))
  tom <- init_tensoromics(est)
  stopifnot(dim(tom$vec_container) == c(8,50),
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
  fn <- file.path(tempdir(), "test_tensoromics.bin")
  save_tensoromics(tom2, fn)
  info <- file.info(fn)
  stopifnot(!is.na(info$size) && info$size > 0)
  cat("Save OK. Datei:", fn, "\n")
}

run_all_tests <- function() {
  cat("=== Testing init & memory ===\n"); test_init_mem()
  cat("=== Testing update       ===\n"); test_update()
  cat("=== Testing save         ===\n"); test_save()
  cat("=== All tests passed! ===\n")
}

# automatisch ausführen, wenn interaktiv
if (interactive()) run_all_tests()

