source("r/tensoromics_functions.R")

# --- Testing for all functions ---
test_array_wrapper <- function(tmpdir = tempdir()) {
  cat("Start Wrapper-Tests...\n")
  fn <- function(name) file.path(tmpdir, name)

  cat("Integer tests...\n")
  cat("Dim 1: ")
  arr1 <- 1:10
  # Debug: check storage mode and length
  stopifnot(is.integer(arr1))
  stopifnot(length(arr1) == 10)
  tox_serialize_int_array(arr1, fn("int1d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int1d.bin")) == arr1))
  cat("okay\n")

  cat("Dim 2: ")
  arr2 <- matrix(1:12, nrow=3, ncol=4)
  stopifnot(is.integer(arr2))
  stopifnot(length(arr2) == 12)
  tox_serialize_int_array(arr2, fn("int2d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int2d.bin")) == arr2))
  cat("okay\n")

  cat("Dim 3: ")
  arr3 <- array(1:24, dim = c(2,3,4))
  tox_serialize_int_array(arr3, fn("int3d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int3d.bin")) == arr3))
  cat("okay\n")

  cat("Dim 4: ")
  arr4 <- array(1:48, dim = c(2,3,4,2))
  tox_serialize_int_array(arr4, fn("int4d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int4d.bin")) == arr4))
  cat("okay\n")

  cat("Dim 5: ")
  arr5 <- array(1:96, dim = c(2,3,4,2,2))
  tox_serialize_int_array(arr5, fn("int5d.bin"))
  stopifnot(all(tox_deserialize_int_array(fn("int5d.bin")) == arr5))
  cat("okay\n")

  # REAL
  cat("Real tests...\n")
  arr1r <- as.numeric(1:10) * 0.5
  tox_serialize_real_array(arr1r, fn("real1d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real1d.bin")) == arr1r))
  cat("okay\n")

  arr2r <- matrix(runif(12), nrow=3, ncol=4)
  tox_serialize_real_array(arr2r, fn("real2d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real2d.bin")) == arr2r))
  cat("okay\n")

  arr3r <- array(runif(24), dim = c(2,3,4))
  tox_serialize_real_array(arr3r, fn("real3d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real3d.bin")) == arr3r))
  cat("okay\n")

  arr4r <- array(runif(48), dim = c(2,3,4,2))
  tox_serialize_real_array(arr4r, fn("real4d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real4d.bin")) == arr4r))
  cat("okay\n")

  arr5r <- array(runif(96), dim = c(2,3,4,2,2))
  tox_serialize_real_array(arr5r, fn("real5d.bin"))
  stopifnot(all(tox_deserialize_real_array(fn("real5d.bin")) == arr5r))
  cat("okay\n")
  cat("Real tests successful!\n")

  # CHAR
  cat("Character tests...\n")
  clen <- 8
  arr1c <- sprintf("%0*d", clen, 1:10)
  tox_serialize_char_array(arr1c, fn("char1d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char1d.bin")) == arr1c))

  arr2c <- matrix(sprintf("%0*d", clen, 1:12), nrow=3, ncol=4)
  tox_serialize_char_array(arr2c, fn("char2d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char2d.bin")) == arr2c))

  arr3c <- array(sprintf("%0*d", clen, 1:24), dim = c(2,3,4))
  tox_serialize_char_array(arr3c, fn("char3d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char3d.bin")) == arr3c))

  arr4c <- array(sprintf("%0*d", clen, 1:48), dim = c(2,3,4,2))
  tox_serialize_char_array(arr4c, fn("char4d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char4d.bin")) == arr4c))

  arr5c <- array(sprintf("%0*d", clen, 1:96), dim = c(2,3,4,2,2))
  tox_serialize_char_array(arr5c, fn("char5d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char5d.bin")) == arr5c))

  # 1D-Array with different char lengths
  arr_ascii1 <- c("A", "G1", "GENE003", "BRCA1", "XYZ", "", "12345678", "SEQ")
  tox_serialize_char_array(arr_ascii1, fn("char_ascii1d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii1d.bin")) == arr_ascii1))

  # 2D-Matrix with different length
  arr_ascii2 <- matrix(c("GENE1", "GENE22", "GENE333", "", "ID", "SEQ9999"), nrow = 2, byrow = TRUE)
  tox_serialize_char_array(arr_ascii2, fn("char_ascii2d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii2d.bin")) == arr_ascii2))

  # 3D-Array with realistic data
  arr_ascii3 <- array(c("TP53", "BRCA1", "MT-ATP6", "CYTB", "ND1", "", "NRAS", "EGFR"), dim = c(2, 2, 2))
  tox_serialize_char_array(arr_ascii3, fn("char_ascii3d.bin"))
  stopifnot(all(tox_deserialize_char_array(fn("char_ascii3d.bin")) == arr_ascii3))

  cat("All ASCII string variation tests passed!\n")
}

# start tests
test_array_wrapper()