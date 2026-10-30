# Generated. Do not edit.

#' The number of threads a parallel loop of this build can use: 1 without OpenMP.
#'
#' A build with OpenMP answers the OpenMP runtime's `omp_get_max_threads()`, which follows
#' the environment variable `OMP_NUM_THREADS` and is otherwise the runtime's own choice,
#' usually the number of processors.
#' A serial build answers 1, as one thread runs every loop there.
#'
#' Generated from the Fortran procedure \code{f42_openmp::max_omp_threads}, whose argument names
#' are the ones an error message reports.
#'
#' @return a integer scalar. The number of threads a parallel loop can use; 1 without OpenMP
#' @export
max_omp_threads <- function() {
    .result <- .Call("max_omp_threads_call")
    .arguments <- c()
    .status <- check_err_code(.result$ierr, .arguments)

    .result$max_threads
}
