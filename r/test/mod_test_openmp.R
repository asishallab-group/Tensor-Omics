# max_omp_threads() reports what the build does with OpenMP: 1 in a serial build, and in an
# OpenMP build the runtime's maximum, which follows OMP_NUM_THREADS.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

test_max_omp_threads_is_a_positive_integer <- function() {
  max_threads <- max_omp_threads()
  assert_true(is.integer(max_threads) && length(max_threads) == 1L, "a single integer")
  assert_true(max_threads >= 1L, "at least one thread")
}

test_max_omp_threads_follows_omp_num_threads <- function() {
  # in a new process: the OpenMP runtime reads OMP_NUM_THREADS once, when it is loaded, and R
  # itself may have loaded it at startup (an OpenBLAS built with libgomp does)
  output <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote('source("r/load_tensor_omics.R"); cat(max_omp_threads())')),
                    env = "OMP_NUM_THREADS=3", stdout = TRUE)
  max_threads <- as.integer(output[length(output)])
  assert_true(max_threads %in% c(1L, 3L), "1 (serial) or OMP_NUM_THREADS=3 (OpenMP)")
}

run_all_tests()
