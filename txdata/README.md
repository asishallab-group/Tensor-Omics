# TensorOmics: Fortran-R Interface for Omics Data Handling

This project provides a high-performance Fortran backend with an R interface for managing and processing tensor-structured omics data. It includes memory-efficient data containers, serialization, and update mechanisms tailored for large-scale biological datasets.

## Table of Contents
- [Project Overview](#project-overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [API Documentation](#api-documentation)
- [Testing](#testing)

---

## Project Overview

### Key Features
- **Memory-Efficient Storage**: Column-major arrays for `vec_container` and `shift_vecs`.
- **Metadata Support**: Structured handling of gene/tissue metadata.
- **Binary Serialization**: Save/load data to `.txdata` files using portable byte streams.
- **Patch Updates**: Incremental updates with index tracking.

### Components
- **Fortran Modules**: Core logic for data structures and computations.
- **R Wrapper**: User-friendly interface for R users.
- **Test Suite**: Validation scripts with data consistency checks.

---

## Prerequisites

- **Fortran Compiler**: `gfortran` (tested with GNU Fortran 10+).
- **R Environment**: R 4.0+ with base packages.
- **Linux/macOS**: Scripts are optimized for Unix-like systems.

---

## Installation

### 1. Compile Fortran Code
```bash
# Navigate to the project directory
cd ~/EasyVectorOmics/tensor-omics/txdata

# Make the compile script executable
chmod +x compile_files.sh

# Compile the Fortran code into a shared library
./compile_files.sh
```

This generates build/libtensoromics.so, which is used by the R wrapper.
Usage
1. Load the R Wrapper
`source("RWrapper.R")`

2. Initialize Data Container

estimates <- create_estimates(
  n_tissues = 10L,
  n_genes = 100L,
  meta_n_rows = 50L,
  meta_max_char = 32L,
  meta_col_types = c(2L, 2L, 3L)  # 2 = real, 3 = char
)

tom <- init_tensoromics(estimates)

3. Update with New Data

patch <- matrix(rnorm(10 * 5), ncol = 5)
tom <- update_tensoromics(tom, patch)

4. Save and Reload Data

save_tensoromics(tom, "data.txdata")

## Read saved data for verification
loaded_data <- read_tensoromics("data.txdata", n_tissues = 10, n_genes = 100)

## API Documentation
### R Functions
Function                                            | Description
create_estimates(...)	                            | Creates metadata estimates for the data container.
init_tensoromics(estimates)	                        | Initializes the tensor container with metadata.
update_tensoromics(tom, patch)	                    | Inserts a data patch and returns updated indices.
save_tensoromics(tom, filename) 	                | Serializes the container to a binary file using portable byte streams.
read_tensoromics(filename, n_tissues, n_genes)	    | Reads binary files into R matrices.
calculate_memory(estimates)	                        | Estimates memory requirements in bytes.

### Fortran Subroutines
1. calculate_memory_requirements: Computes memory usage.
2. init: Allocates and initializes data structures.
3. update: Handles patch insertion and index tracking.
4. save: Writes raw binary data using stream I/O (no record markers).

# Testing
## Run the built-in test suite:
### In R:
run_all_tests()

Expected Output:

=== Testing init & memory ===
Init OK. Memory: 6984 bytes
=== Testing update       ===
Update OK. Indices: 1 2 3 4 5 
=== Testing save         ===
Save OK. File: data/test_tensoromics.txdata 
Data consistency verified!
=== All tests passed! ===
