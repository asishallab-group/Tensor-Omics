# max_omp_threads() reports what the build does with OpenMP: 1 in a serial build, and in an
# OpenMP build the runtime's maximum, which follows OMP_NUM_THREADS.
#
# What the build was MEANT to be comes from TOX_EXPECT_OPENMP, which run_all_tests.sh sets from its
# --openmp option (1 with it, 0 without). Run by hand without the variable, the test accepts either
# answer, as it cannot tell which build it runs against.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_THREADS_REQUESTED <- 3L

# TOX_EXPECT_OPENMP -> the answers max_omp_threads() may give under OMP_NUM_THREADS=3; unset, either
ALLOWED_ANSWERS <- list("1" = N_THREADS_REQUESTED, "0" = 1L)

test_max_omp_threads_is_a_positive_integer <- function() {
  max_threads <- max_omp_threads()
  assert_true(is.integer(max_threads) && length(max_threads) == 1L, "a single integer")
  assert_true(max_threads >= 1L, "at least one thread")
}

test_max_omp_threads_follows_omp_num_threads <- function() {
  expect_openmp <- Sys.getenv("TOX_EXPECT_OPENMP")
  allowed <- if (expect_openmp == "") c(1L, N_THREADS_REQUESTED) else ALLOWED_ANSWERS[[expect_openmp]]
  assert_true(!is.null(allowed), sprintf("TOX_EXPECT_OPENMP is 1, 0 or unset, not '%s'", expect_openmp))
  # in a new process: the OpenMP runtime reads OMP_NUM_THREADS once, when it is loaded, and R
  # itself may have loaded it at startup (an OpenBLAS built with libgomp does)
  output <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote('source("r/load_tensor_omics.R"); cat(max_omp_threads())')),
                    env = paste0("OMP_NUM_THREADS=", N_THREADS_REQUESTED), stdout = TRUE)
  max_threads <- as.integer(output[length(output)])
  assert_true(max_threads %in% allowed,
              sprintf("%s with TOX_EXPECT_OPENMP='%s' and OMP_NUM_THREADS=%d, not %s",
                      paste(allowed, collapse = " or "), expect_openmp, N_THREADS_REQUESTED, max_threads))
}

run_all_tests()
