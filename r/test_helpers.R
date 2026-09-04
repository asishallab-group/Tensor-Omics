ANSI_START <- "\033["
FG256_START <- paste0(ANSI_START, "38;5;")

COLORS <- list(
  green = 154,
  copper = 214,
  dark_copper = 208,
  red = 196,
  light_gray = 252,
  yellow = 226,
  cream = 255,
  error = 222
)

fg256 <- function(color_name = "") {
  if (!is.null(COLORS[[color_name]])) {
    paste0(FG256_START, COLORS[[color_name]], "m")
  } else {
    paste0(ANSI_START, "0m")
  }
}

ccat <- function(...) {
  colored <- paste0(...)
  parts <- strsplit(colored, "@")[[1]]

  for (substr in parts) {
    color_name <- strsplit(substr, "\\.", fixed = FALSE)[[1]][1]
    pattern <- paste0("@", color_name, "\\.")
    replacement <- fg256(color_name)
    colored <- sub(pattern, replacement, colored, fixed = FALSE)
  }

  cat(colored, fg256(), sep="")
}

run_all_tests <- function(env = parent.frame(), test_only = TRUE) {
  # Discover candidate names
  env <- as.environment(env)
  if (test_only) {
    test_names <- ls(env, pattern = "^test_")
  } else {
    test_names <- ls(env)
  }

  passed  <- 0
  failed  <- 0
  skipped <- 0

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_name <- sprintf("@light_gray.%s@light_gray.", sub("^--file=", "", file_arg))

  ccat("@cream.Running R tests of '", script_name, "@cream.'...\n")

  for (name in test_names) {
    obj <- get(name, envir = env)
    if (!is.function(obj)) next  # skip non-functions

    test_func <- obj

    tryCatch(
      {

        test_func()
        ccat(sprintf("@green.✓ @copper.%s @green.passed@cream..\n", name))
        passed <- passed + 1
      },
      error = function(e) {
        msg <- conditionMessage(e)

        # # Same skip logic as Python version
        # if (grepl("Note:", msg) || grepl("acceptable", msg, ignore.case = TRUE)) {
        #   ccat(sprintf("~ %s skipped (expected behavior): %s\n", name, msg))
        #   skipped <<- skipped + 1
        # } else {
          ccat(sprintf("@red.✗ @dark_copper.%s @red.FAILED@cream.: @error.%s\n", name, msg))
          failed <<- failed + 1
        # }
      }
    )
  }

  ccat("@cream.\nSummary: @green.", passed, " passed@cream., @red.", failed, " failed@cream., @yellow.", skipped, " skipped\n")

  if (failed) quit(status = 1)

  ccat("@cream.All tests in '", script_name, "@cream.' passed successfully.\n")
}

assert_true <- function(expr, msg = "Assertion failed") {
  if (!isTRUE(expr)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_false <- function(expr, msg = "Assertion failed") {
  if (isTRUE(expr)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

.tox_code_label <- function(code) {
  code <- as.integer(code)
  for (name in ls(globalenv(), pattern = "^ERR_")) {
    value <- get(name, envir = globalenv())
    if (is.numeric(value) && length(value) == 1L && as.integer(value) == code)
      return(sprintf("%s (%d)", name, code))
  }
  sprintf("code %d", code)
}

# `code` is a bare tox_errors code (ERR_INVALID_INPUT), not a packed ierr: check_err_code
# strips the argument position before it signals. Omit it to demand an error carrying no tox
# code -- a plain R error, or one of the binding's own tox_type_error / tox_shape_error
# conditions, which are tox_error conditions with code = NA by design.
assert_error <- function(expr, msg = "Expected an error", code = NULL) {
  if (!is.null(code) && as.integer(code) >= ARG_POS_FACTOR)
    stop(sprintf("assert_error takes a bare error code, not a packed ierr: %d", as.integer(code)),
         call. = FALSE)

  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(msg, ": nothing was raised", call. = FALSE)

  actual <- if (inherits(err, "tox_error")) err$code else NULL
  actual <- if (is.null(actual) || length(actual) != 1L) NA_integer_ else as.integer(actual)

  if (is.null(code)) {
    if (!is.na(actual))
      stop(msg, ": expected an error that is not a tensor-omics error, but got ",
           .tox_code_label(actual), ": ", conditionMessage(err), call. = FALSE)
  } else if (is.na(actual) || actual != as.integer(code)) {
    got <- if (is.na(actual)) sprintf("<%s>", paste(class(err), collapse = ", "))
           else .tox_code_label(actual)
    stop(msg, ": expected ", .tox_code_label(code), ", but got ", got, ": ",
         conditionMessage(err), call. = FALSE)
  }
  invisible(TRUE)
}

assert_equal_int <- function(x, y, msg = "Integer mismatch") {
  if (!identical(x, y)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

assert_equal_numeric <- function(x, y, tol = 1e-12, msg = "Real mismatch") {
  if (length(x) != length(y)) stop(msg, call. = FALSE)
  for (i in seq_len(length(x))) {
    if (!(is.finite(x[i]) && is.finite(y[i]) && abs(x[i] - y[i]) <= tol)) stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}