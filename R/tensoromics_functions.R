# tensoromics_functions.R

# --- Library Loading ---
lib_path <- file.path("build", "libtensor-omics.so")
if (!file.exists(lib_path)) {
  stop(paste("Shared library not found at:", lib_path,
             "\nPlease build the Fortran project using './build.sh' first."))
}
dyn.load(lib_path)

# This constant must match the Fortran parameter
MAX_FIELD_LEN <- 512

# --- Error Handling ---
check_fortran_errors <- function(status_code, context) {
  msg <- switch(
    as.character(status_code),
    "0" = NULL,
    "1" = "No data to process.",
    "2" = "Column index out of bounds.",
    "3" = "At least one data conversion error occurred.",
    "4" = "The specified file is empty.",
    "5" = "File has a header but no data rows.",
    "10" = "File not found or cannot be opened.",
    paste("An unknown Fortran error occurred (code:", status_code, ")")
  )
  if (!is.null(msg)) {
    stop(paste("Error in", context, ":", msg), call. = FALSE)
  }
}

# --- Utilities ---
.string_to_ascii <- function(s) {
  if (nchar(s) == 0) return(integer(0))
  as.integer(charToRaw(s))
}

# --- User-Facing API Functions ---

#' Get CSV dimensions without reading data.
#' @export
tox.read_csv_dimensions <- function(filepath, has_header = TRUE, delimiter = ",") {
  stopifnot(file.exists(filepath), is.character(filepath),
            is.logical(has_header), nchar(delimiter) == 1)
            
  fname_ascii <- .string_to_ascii(filepath)
  res <- .Fortran("get_csv_dims_r",
                   filename_ascii = as.integer(fname_ascii),
                   fn_len = as.integer(length(fname_ascii)),
                   has_header = as.logical(has_header),
                   delimiter_ascii = as.integer(charToRaw(delimiter)),
                   num_rows = integer(1),
                   num_cols = integer(1),
                   status = integer(1))
  
  check_fortran_errors(res$status, "tox.read_csv_dimensions")
  list(num_rows = res$num_rows, num_cols = res$num_cols)
}

#' Reads a CSV file, returning raw ASCII arrays for header and data.
#' @export
tox.read_csv_raw <- function(filepath, has_header = TRUE, delimiter = ",") {
  dims <- tox.read_csv_dimensions(filepath, has_header, delimiter)
  fname_ascii <- .string_to_ascii(filepath)
  
  header_out <- if (has_header) array(0L, dim = c(MAX_FIELD_LEN, dims$num_cols)) else array(0L, dim = c(0,0))
  data_out <- array(0L, dim = c(MAX_FIELD_LEN, dims$num_rows, dims$num_cols))
  
  res <- .Fortran("read_csv_to_strings_r",
                   filename_ascii = as.integer(fname_ascii),
                   fn_len = as.integer(length(fname_ascii)),
                   has_header = as.logical(has_header),
                   delimiter_ascii = as.integer(charToRaw(delimiter)),
                   num_rows = as.integer(dims$num_rows),
                   num_cols = as.integer(dims$num_cols),
                   header_out_ascii = header_out,
                   data_out_ascii = data_out,
                   status = integer(1))
                   
  check_fortran_errors(res$status, "tox.read_csv_raw")
  
  list(
    header = if (has_header) res$header_out_ascii else NULL,
    data = res$data_out_ascii,
    dims = dims
  )
}

#' Converts a raw ASCII array from Fortran into an R character matrix or vector.
#' @export
tox.convert_ascii_to_strings <- function(ascii_array) {
  if (is.null(ascii_array) || length(ascii_array) == 0) return(NULL)
  
  dims <- dim(ascii_array)
  is_3d <- length(dims) == 3
  
  num_rows <- if (is_3d) dims[2] else 1
  num_cols <- if (is_3d) dims[3] else dims[2]
  mat <- matrix("", nrow = num_rows, ncol = num_cols)
  
  for (j in 1:num_cols) {
    for (i in 1:num_rows) {
      raw_vec <- if (is_3d) ascii_array[, i, j] else ascii_array[, j]
      null_pos <- which(raw_vec == 0)
      end_pos <- if (length(null_pos) > 0) null_pos[1] - 1 else length(raw_vec)
      if (end_pos > 0) {
        mat[i, j] <- rawToChar(as.raw(raw_vec[1:end_pos]))
      }
    }
  }
  return(if (num_rows == 1) as.vector(mat) else mat)
}

#' Generic internal function to convert raw ASCII data to a numeric type.
.read_typed_columns <- function(raw_data, columns, type_name, fortran_func_name) {
  stopifnot(is.array(raw_data), is.numeric(columns), all(columns > 0))
  
  dims <- dim(raw_data)
  num_rows <- dims[2]
  num_cols_in <- dims[3]
  num_cols_to_read <- length(columns)
  
  output_matrix <- switch(
    type_name,
    "integer" = matrix(0L, nrow = num_rows, ncol = num_cols_to_read),
    "double"  = matrix(0.0, nrow = num_rows, ncol = num_cols_to_read),
    "logical" = matrix(FALSE, nrow = num_rows, ncol = num_cols_to_read),
    "complex" = matrix(complex(real=0, imaginary=0), nrow = num_rows, ncol = num_cols_to_read)
  )
  
  res <- .Fortran(fortran_func_name,
                  data_in_ascii = raw_data,
                  num_rows = as.integer(num_rows),
                  num_cols_in = as.integer(num_cols_in),
                  cols_to_read = as.integer(columns),
                  num_cols_to_read = as.integer(num_cols_to_read),
                  output_data = output_matrix,
                  status = integer(1),
                  # Rename the output argument in the wrapper to a consistent name
                  PACKAGE = "tensoromics"
                 )
  
  check_fortran_errors(res$status, fortran_func_name)
  res$output_data
}

# The Fortran wrappers need to be updated to have a consistent output argument name, e.g., "output_data"
# For example: SUBROUTINE read_integer_columns_r(..., int_data_out, status) ->
# SUBROUTINE read_integer_columns_r(..., output_data, status)
# For now, we call them individually.

#' @export
tox.read_integer_columns <- function(raw_data, columns) {
  res <- .Fortran("read_integer_columns_r", data_in_ascii=raw_data, num_rows=as.integer(dim(raw_data)[2]), num_cols_in=as.integer(dim(raw_data)[3]), cols_to_read=as.integer(columns), num_cols_to_read=as.integer(length(columns)), int_data_out=matrix(0L, nrow=dim(raw_data)[2], ncol=length(columns)), status=integer(1))
  check_fortran_errors(res$status, "read_integer_columns_r")
  res$int_data_out
}

#' @export
tox.read_real_columns <- function(raw_data, columns) {
  res <- .Fortran("read_real_columns_r", data_in_ascii=raw_data, num_rows=as.integer(dim(raw_data)[2]), num_cols_in=as.integer(dim(raw_data)[3]), cols_to_read=as.integer(columns), num_cols_to_read=as.integer(length(columns)), real_data_out=matrix(0.0, nrow=dim(raw_data)[2], ncol=length(columns)), status=integer(1))
  check_fortran_errors(res$status, "read_real_columns_r")
  res$real_data_out
}

#' @export
tox.read_logical_columns <- function(raw_data, columns) {
  res <- .Fortran("read_logical_columns_r", data_in_ascii=raw_data, num_rows=as.integer(dim(raw_data)[2]), num_cols_in=as.integer(dim(raw_data)[3]), cols_to_read=as.integer(columns), num_cols_to_read=as.integer(length(columns)), logical_data_out=matrix(FALSE, nrow=dim(raw_data)[2], ncol=length(columns)), status=integer(1))
  check_fortran_errors(res$status, "read_logical_columns_r")
  res$logical_data_out
}

#' @export
tox.read_complex_columns <- function(raw_data, columns) {
  res <- .Fortran("read_complex_columns_r", data_in_ascii=raw_data, num_rows=as.integer(dim(raw_data)[2]), num_cols_in=as.integer(dim(raw_data)[3]), cols_to_read=as.integer(columns), num_cols_to_read=as.integer(length(columns)), complex_data_out=matrix(complex(1), nrow=dim(raw_data)[2], ncol=length(columns)), status=integer(1))
  check_fortran_errors(res$status, "read_complex_columns_r")
  res$complex_data_out
}

#' @export
tox.read_character_columns <- function(raw_data, columns) {
  num_rows <- dim(raw_data)[2]
  # Select columns and convert to string matrix
  tox.convert_ascii_to_strings(raw_data[,,columns, drop=FALSE])
}