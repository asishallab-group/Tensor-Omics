# Comprehensive R test suite for LOESS interface functions
# Validation of the R wrapper for plain and robust LOESS.

source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

# 1. Test Workspace Calculation
test_workspace_calculation <- function() {
  
  # Typical parameters: d=1 (univariate), nvmax=100, setlf=TRUE
  ws <- tox_loess_required_workspace(d=1, nvmax=100, setlf=TRUE)
  
  assert_true(is.list(ws))
  assert_true(ws$liv > 0)
  assert_true(ws$lv > 0)
  
}

# 2. Test Plain LOESS Functionality
test_loess_plain_functionality <- function() {
  
  n <- 20
  x <- seq(1, 10, length.out = n)
  # Generate a linear trend with minor noise
  y <- 2.0 * x + rnorm(n, mean = 0, sd = 0.1)
  w <- rep(1.0, n)
  z <- x
  
  # Get required workspace sizes
  ws <- tox_loess_required_workspace(d=1, nvmax=n, setlf=FALSE)
  
  # Initialize workspace arrays
  iv <- integer(ws$liv)
  wv <- numeric(ws$lv)
  
  yhat <- loess_fit_plain(
    x=x, y=y, w=w, z=z, 
    span=0.5, degree=1, nvmax=n, 
    infl=FALSE, setlf=FALSE, 
    iv=iv, wv=wv
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
  ws <- tox_loess_required_workspace(d=1, nvmax=n, setlf=FALSE)
  
  # Workspace and additional robust arrays
  iv <- integer(ws$liv)
  wv <- numeric(ws$lv)
  rw <- numeric(n)
  ww <- numeric(n)
  res <- numeric(n)
  pi <- integer(n)
  
  yhat <- loess_fit_robust(
    x=x, y=y, w=w, z=z, 
    span=0.5, degree=1, nvmax=n, 
    infl=FALSE, setlf=FALSE, n_iters=4,
    iv=iv, wv=wv, rw=rw, ww=ww, res=res, pi=pi
  )
  
  assert_true(length(yhat) == n)
  # If robustness works, the outlier at index 6 should be ignored
  assert_true(yhat[6] < 50.0)
  
}

# 4. Test High-level Wrapper (tox_loess)
test_tox_loess_wrapper <- function() {
  
  n <- 30
  x <- seq(0, 2*pi, length.out = n)
  y <- sin(x)
  
  # Test Plain mode (mode=0)
  yhat_plain <- tox_loess(x, y, span=0.4, degree=1, mode=0)
  assert_true(length(yhat_plain) == n)
  
  # Test Robust mode (mode=1)
  yhat_robust <- tox_loess(x, y, span=0.4, degree=1, mode=1, n_iters=2)
  assert_true(length(yhat_robust) == n)
  
  # Verify that results are not identical due to re-weighting
  assert_true(!all(yhat_plain == yhat_robust))
  
}

# 5. Test Error Handling
test_invalid_inputs <- function() {
  
  x <- 1:10
  y <- 1:5 # Mismatched length
  
  error_caught <- FALSE
  tryCatch({
    tox_loess(x, y)
  }, error = function(e) {
    error_caught <<- TRUE
  })
  
  assert_true(error_caught)
}

run_all_tests()