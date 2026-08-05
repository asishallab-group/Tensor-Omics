# Comprehensive R test suite for LOESS binding functions
# Validation of the R wrapper for plain and robust LOESS.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# 1. Test Workspace Calculation
test_workspace_calculation <- function() {
  
  # Typical parameters: univariate, 100 neighbours, factorization saved
  ws <- tox_loess_required_workspace(n_dim=1, max_neighborhood_size=100, save_factorization=TRUE)
  
  assert_true(is.list(ws))
  assert_true(ws$int_workspace_size > 0)
  assert_true(ws$real_workspace_size > 0)
  
}

# 2. Test Plain LOESS Functionality
test_loess_plain_functionality <- function() {
  
  n <- 20
  x <- seq(1, 10, length.out = n)
  # Generate a linear trend with minor noise
  y <- 2.0 * x + rnorm(n, mean = 0, sd = 0.1)
  w <- rep(1.0, n)
  z <- x
  
  yhat <- loess_fit_plain(
    x=x, y=y, weights=w, eval_points=matrix(z, ncol=1),
    span=0.5, degree=1, max_neighborhood_size=n,
    compute_influence=FALSE, save_factorization=FALSE
  )
  
  assert_true(length(yhat) == n)
  assert_true(all(!is.na(yhat)))
  
}

# 3. Test Robust LOESS Functionality (Outlier Suppression)
test_loess_robust_functionality <- function() {
  
  n <- 20
  x <- seq(1, 10, length.out = n)
  y <- 3.0 * x
  y[6] <- 100.0  # Introduce an aggressive outlier (index 6 in R is 5 in 0-based)
  
  w <- rep(1.0, n)
  z <- x

  # The robust scratch arrays (robust/combined weights, residuals, permutation) are now
  # allocated by the wrapper, so the friendly binding no longer takes them.
  yhat <- loess_fit_robust(
    x=x, y=y, weights=w, eval_points=matrix(z, ncol=1),
    span=0.5, degree=1, max_neighborhood_size=n,
    compute_influence=FALSE, save_factorization=FALSE, n_iters=4
  )
  
  assert_true(length(yhat) == n)
  # If robustness works, the outlier at index 6 should be ignored
  assert_true(yhat[6] < 50.0)
  
}

# 4. The plain and robust fits are separate entry points
test_plain_and_robust_are_separate_entry_points <- function() {

  n <- 30
  x <- seq(0, 2*pi, length.out = n)
  y <- sin(x)
  y[8] <- y[8] + 5.0  # an outlier for the robust iterations to down-weight

  # the self-allocating entry points still take what the fit is *of*: uniform weights,
  # evaluated at the training x, is the common case
  w <- rep(1.0, n)
  z <- matrix(x, ncol=1)

  yhat_plain <- loess_fit_plain(x=x, y=y, weights=w, eval_points=z,
                                span=0.4, degree=1, max_neighborhood_size=n)
  assert_true(length(yhat_plain) == n)

  yhat_robust <- loess_fit_robust(x=x, y=y, weights=w, eval_points=z,
                                  span=0.4, degree=1, max_neighborhood_size=n, n_iters=2)
  assert_true(length(yhat_robust) == n)

  # the robust fit suppresses the outlier, so the two cannot agree
  assert_true(!all(yhat_plain == yhat_robust))

}

# 5. Test Error Handling
test_invalid_inputs <- function() {

  n <- 10
  x <- seq(1, 10, length.out = n)
  y <- x
  w <- rep(1.0, n)
  z <- matrix(x, ncol=1)

  error_caught <- FALSE
  tryCatch({
    # a span outside (0, 1] is rejected by the generated wrapper
    loess_fit_plain(x=x, y=y, weights=w, eval_points=z,
                    span=2.0, degree=1, max_neighborhood_size=n)
  }, error = function(e) {
    error_caught <<- TRUE
  })

  assert_true(error_caught)
}

run_all_tests()