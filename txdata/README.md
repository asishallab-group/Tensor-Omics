```markdown
# TensorOmics: Fortran-R Interface for Omics Data Handling

This project provides a high-performance Fortran backend with an R interface for managing and processing tensor-structured omics data. It includes memory-efficient data containers, serialization, and update mechanisms tailored for large-scale biological datasets.

## Table of Contents
- [Project Overview](#project-overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [API Documentation](#api-documentation)
- [Testing](#testing)
- [License](#license)

---

## Project Overview

### Key Features
- **Memory-Efficient Storage**: Column-major arrays for `vec_container` and `shift_vecs`.
- **Metadata Support**: Structured handling of gene/tissue metadata.
- **Binary Serialization**: Save/load data to `.txdata` files.
- **Patch Updates**: Incremental updates with index tracking.

### Components
- **Fortran Modules**: Core logic for data structures and computations.
- **R Wrapper**: User-friendly interface for R users.
- **Test Suite**: Validation scripts for functionality.

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

This generates `libtensoromics.so`, which is used by the R wrapper.

---

## Usage

### 1. Load the R Wrapper
```r
source("RWrapper.R")
```

### 2. Initialize Data Container
```r
estimates <- create_estimates(
  n_tissues = 10L,
  n_genes = 100L,
  meta_n_rows = 50L,
  meta_max_char = 32L,
  meta_col_types = c(2L, 2L, 3L)  # 2 = real, 3 = char

tom <- init_tensoromics(estimates)
```

### 3. Update with New Data
```r
patch <- matrix(rnorm(10 * 5), ncol = 5)
tom <- update_tensoromics(tom, patch)
```

### 4. Save to File
```r
save_tensoromics(tom, "data.txdata")
```

---

## API Documentation

### R Functions

| Function                     | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| `create_estimates(...)`      | Creates metadata estimates for the data container.                          |
| `init_tensoromics(estimates)`| Initializes the tensor container with metadata.                             |
| `update_tensoromics(tom, patch)` | Inserts a data patch and returns updated indices.                       |
| `save_tensoromics(tom, filename)` | Serializes the container to a binary file.                              |
| `calculate_memory(estimates)`| Estimates memory requirements in bytes.                                     |

### Fortran Subroutines
- `calculate_memory_requirements`: Computes memory usage.
- `init`: Allocates and initializes data structures.
- `update`: Handles patch insertion and index tracking.
- `save`: Writes data to binary format.

---

## Testing

Run the built-in test suite:
```r
# In R:
run_all_tests()
```

**Expected Output**:
```
=== Testing init & memory ===
Init OK. Memory: [bytes] 
=== Testing update       ===
Update OK. Indices: 1 2 3 4 5 
=== Testing save         ===
Save OK. Datei: /tmp/test_tensoromics.bin 
=== All tests passed! ===
```

! CURRENTLY THE .Fortran() CALL DOES NOT WORK !
