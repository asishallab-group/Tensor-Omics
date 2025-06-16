# On Windows, use "libF42CsvReader.dll"
# On macOS/Linux, use "libF42CsvReader.so"
shared_lib_file <- "build/libF42CsvReader.so"
if (.Platform$OS.type == "windows") {
  shared_lib_file <- "build/libF42CsvReader.dll"
}

# Load the shared library
dyn.load(shared_lib_file)

# --- Define R wrapper functions that use the modern .Call() interface ---

read_csv <- function(filename, has_header = TRUE, delimiter = ",") {
  # Call the C-bound Fortran routine. Note the "_C" suffix.
  result <- .Call("read_csv_file_C",
                  as.character(filename),
                  as.logical(has_header),
                  as.character(delimiter))
  return(result)
}

serialize_data <- function(filename) {
  result <- .Call("serialize_C", as.character(filename))
  return(result)
}

deserialize_data <- function(filename) {
  result <- .Call("deserialize_C", as.character(filename))
  return(result)
}


# --- Run Tests ---

cat("--> Testing C binding via R's .Call() interface\n")

io_status <- read_csv("test_data.csv", has_header = TRUE, delimiter = ",")
if (io_status != 0) {
  stop(paste("FAILURE: read_csv_file_C failed with status", io_status))
}
cat("SUCCESS: read_csv_file_C passed.\n")

io_status <- serialize_data("test_data.bin")
if (io_status != 0) {
  stop(paste("FAILURE: serialize_C failed with status", io_status))
}
cat("SUCCESS: serialize_C passed.\n")

io_status <- deserialize_data("test_data.bin")
if (io_status != 0) {
  stop(paste("FAILURE: deserialize_C failed with status", io_status))
}
cat("SUCCESS: deserialize_C passed.\n")