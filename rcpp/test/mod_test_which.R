source("rcpp/tensoromics_functions.R")

cat("=== Testing tox_which R wrapper ===\n")

# Test 1: Basic integer mask
mask1 <- c(0L, 1L, 0L, 1L, 1L)
res1 <- tox_which(mask1)
stopifnot(identical(as.integer(res1), c(2L, 4L, 5L)))
cat("test_which_basic_integer passed ✓\n")

# Test 2: Logical mask input
mask2 <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
res2 <- tox_which(mask2)
stopifnot(identical(as.integer(res2), c(2L, 4L)))
cat("test_which_logical_input passed ✓\n")

# Test 3: Respect m_max cap
mask3 <- c(1L, 0L, 1L, 1L, 0L)
res3 <- tox_which(mask3, m_max = 2L)
stopifnot(length(res3) == 2L)
stopifnot(identical(as.integer(res3), c(1L, 3L)))
cat("test_which_m_max_cap passed ✓\n")

# Test 4: No matches
mask4 <- c(0L, 0L, 0L)
res4 <- tox_which(mask4)
stopifnot(length(res4) == 0L)
cat("test_which_no_matches passed ✓\n")

cat("All tox_which tests passed.\n")
