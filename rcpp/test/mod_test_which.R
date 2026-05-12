source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

test_which <- function() {
    # Test 1: Basic integer mask
    mask1 <- c(0L, 1L, 0L, 1L, 1L)
    res1 <- tox_which(mask1)
    assert_true(identical(as.integer(res1), c(2L, 4L, 5L)))

    # Test 2: Logical mask input
    mask2 <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
    res2 <- tox_which(mask2)
    assert_true(identical(as.integer(res2), c(2L, 4L)))

    # Test 3: Respect m_max cap
    mask3 <- c(1L, 0L, 1L, 1L, 0L)
    res3 <- tox_which(mask3, m_max = 2L)
    assert_true(length(res3) == 2L)
    assert_true(identical(as.integer(res3), c(1L, 3L)))

    # Test 4: No matches
    mask4 <- c(0L, 0L, 0L)
    res4 <- tox_which(mask4)
    assert_true(length(res4) == 0L)
}

run_all_tests()
