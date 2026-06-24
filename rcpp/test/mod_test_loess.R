# Comprehensive R test suite for LOESS interface functions
# Validation of the R wrapper for plain and robust LOESS.

source("rcpp/tensoromics_functions.R")

# 1. Test Workspace Calculation
test_workspace_calculation <- function() {
  cat("Testing tox_loess_required_workspace...\n")
  
  # Typical parameters: d=1 (univariate), nvmax=100, setlf=TRUE
  ws <- tox_loess_required_workspace(d=1, nvmax=100, setlf=TRUE)
  
  stopifnot(is.list(ws))
  stopifnot(ws$liv > 0)
  stopifnot(ws$lv > 0)
  
  cat(sprintf("  Workspace calculation passed (liv=%d, lv=%d)\n", ws$liv, ws$lv))
}

# 2. Test Plain LOESS Functionality
test_loess_plain_functionality <- function() {
  cat("Testing loess_fit_plain...\n")
  
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
  
  stopifnot(length(yhat) == n)
  stopifnot(all(!is.na(yhat)))
  
  cat("  loess_fit_plain passed.\n")
}

# 3. Test Robust LOESS Functionality (Outlier Suppression)
test_loess_robust_functionality <- function() {
  cat("Testing loess_fit_robust...\n")
  
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
  
  stopifnot(length(yhat) == n)
  # If robustness works, the outlier at index 6 should be ignored
  stopifnot(yhat[6] < 50.0)
  
  cat(sprintf("  loess_fit_robust passed (yhat[6]=%.2f vs y[6]=100)\n", yhat[6]))
}

# 4. Test High-level Wrapper (tox_loess)
test_tox_loess_wrapper <- function() {
  cat("Testing tox_loess (High-level wrapper)...\n")
  
  n <- 30
  x <- seq(0, 2*pi, length.out = n)
  y <- sin(x)
  
  # Test Plain mode (mode=0)
  yhat_plain <- tox_loess(x, y, span=0.4, degree=1, mode=0)
  stopifnot(length(yhat_plain) == n)
  
  # Test Robust mode (mode=1)
  yhat_robust <- tox_loess(x, y, span=0.4, degree=1, mode=1, n_iters=2)
  stopifnot(length(yhat_robust) == n)
  
  # Verify that results are not identical due to re-weighting
  stopifnot(!all(yhat_plain == yhat_robust))
  
  cat("  tox_loess wrapper passed.\n")
}

# 5. Test Error Handling
test_invalid_inputs <- function() {
  cat("Testing error handling...\n")
  
  x <- 1:10
  y <- 1:5 # Mismatched length
  
  error_caught <- FALSE
  tryCatch({
    tox_loess(x, y)
  }, error = function(e) {
    error_caught <<- TRUE
  })
  
  stopifnot(error_caught)
  cat("  Error handling passed.\n")
}

# Run all tests
cat("=================================================\n")
cat("      LOESS INTERFACE FULL R TESTS\n")
cat("=================================================\n\n")

test_workspace_calculation()
test_loess_plain_functionality()
test_loess_robust_functionality()
test_tox_loess_wrapper()
test_invalid_inputs()

cat("\n=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all LOESS R interface tests passed successfully!\n")