


# Tensor Omics

## **Table of Contents**
- **[Introduction](#introduction)**
- **[Key Features](#key-features)**
- **[Project Structure](#project-structure)**
- **[Compilation](#compilation)**
  - [Compiler Requirements](#compiler-requirements)
  - [Local Compilation (No Docker)](#local-compilation-no-docker)
  - [Compilation Using Docker](#compilation-using-docker)
- **[Installation](#installation)**
  - [Linux](#linux)
  - [Windows-11](#windows-11)
  - [macOS](#macos)
- **[Installation for R and Python](#installation-for-r-and-python)**
  - [R Integration](#r-integration)
  - [Python Integration](#python-integration)
- **[Usage](#usage)**
- **[Testing](#testing)**
- **[Configuration](#configuration)**
- **[License](#license)**

---

# Introduction

Tensor Omics is a high-performance framework for explainable, geometry-based analysis of multimodal omics and related high-dimensional datasets. Instead of relying on black-box models, it treats expression profiles, clinical measures, or socioeconomic indicators as vectors in semantically meaningful spaces (e.g. tissues, disease stages, conditions). By measuring distances, angles, projections, and trajectories in these spaces, Tensor Omics enables direct comparison of activity across genes, paralogs, sexes, species, or patient groups. This geometric approach makes complex multivariate patterns interpretable and reproducible while remaining robust to sparsity and noise. Designed for distributed high-performance computing, Tensor Omics is implemented in Fortran and C with OpenMP parallelisation, SIMD optimisation, and Fortran Coarrays, making the algorithms embarrassingly parallel and suitable for federated datasets where privacy and efficiency are critical. Scientific use cases include: detecting disease biomarkers and subtype-specific trajectories in medical data; quantifying divergence and neofunctionalization of gene duplicates in plant and animal transcriptomes; and reconstructing global gender-equality trajectories from socioeconomic indicators. Across these domains, Tensor Omics provides a unified, geometry-driven methodology for discovering explanatory patterns in heterogeneous, high-dimensional data.


## Key Features

* **Geometry-based analysis**: distances, angles, projections, and trajectory shifts are used as primary primitives.
* **Explainable outputs**: results are interpretable in terms of vector geometry rather than opaque model coefficients.
* **Multi-modal integration**: unifies transcriptomics, proteomics, metabolomics, clinical, or socioeconomic data within one framework.
* **Parallel and federated**: implemented in Fortran/C with OpenMP, SIMD, and coarrays for efficient large-scale computation on distributed datasets.
* **Robust to sparsity and noise**: percentile-based empirical thresholds and local geometric measures enable stability.
* **Broad applications**: demonstrated on medical biomarker discovery, gene duplication outcomes, developmental trajectories, and socioeconomic indicators.


*Tensor Omics shows that geometry, when treated not as preprocessing but as the central instrument of analysis, can open entirely new ways to read complex biological and social data — simple, transparent, and surprisingly powerful.*

---

# Project Structure

This repository contains the source code, methods, snippets and tests for the **Tensor Omics (TOX)** project.

## Folder Overview

```
/build
  └── ...       # Compiled Fortran binaries and intermediate build files

/doc
  └── ...       # Documentation generated automatically using FORD (Fortran documentation tool)

/misc
  └── ...       # Tensor Omics documentation, coding guides and helper dockerfile to compile the project

/python
  └── ...       # Python scripts that execute pipeline logic and invoke subroutines

/r
  └── ...       # R scripts that execute pipeline logic and invoke subroutines

/snippets
└── ...         # Code templates or reusable short logic blocks 

/src
  └── ...       # Fortran backend

/test
  └── ...       # Fortran testing

/helper
  └── helper_c_wrapper.py       #  helper script to generate c wrappers subroutines

build.sh        # Compile and generate shared libraries
ford.yml        # Generates documentation
fpm.toml        # Defines compilation options
test_runner.sh  # Compile and generate unit test

```

## Notes

* **`/build`** is used to store shared libraries, compiled and binary files resulting from Fortran compilation. It keeps the repo clean by separating source and compiled code.
* **`/doc`** contains the auto-generated documentation, which is built using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) from annotated Fortran source files.
* **`/misc`** contains the team's coding guidelines at [Fortran_Coding_Guides.pdf](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/main/misc/Fortran_Coding_Guides.pdf?ref_type=heads), the detailed description of Tensor Omics at [Tensor_Omics_Methods.pdf](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/main/misc/Tensor_Omics_Methods.pdf?ref_type=heads), and a [Dockerfile](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/main/misc/gfortran.docker?ref_type=heads) to compile the project without needing to install anything except Docker.

* **`/python`** includes python scripts that coordinate analysis workflows
* **`/r`** includes R scripts that coordinate analysis workflows
* **`/snippets/`** includes frequently used or testable units of logic reused across development stages. 
  - Snippets should be easy to create and use. The goal is to give the user access to the subroutine names along with their respective arguments, and nothing more. Example:
  ```
    {
        "Call to subroutine_name": {
          "prefix": "tox|f42:subroutine_name",
          "body": [
            "! Brief explanation on what the subroutine does",
            "call subroutine_name(arg1, arg2, arg3, arg4)"
          ],
          "description": "Insert a call to subroutine_name with brief explanation"
        }
    }
  ```
* **`/src`** contains performance-critical Fortran code. These are compiled during the build process.
  - All `.f90` files should include `precompiler_constants.f90`
  - Subroutines that do not perform `input/output` operations or memory allocations must be declared as `pure`.
* **`/test`** contains the unit tests for the Fortran subroutines.

  * The file `asserts.f90` must exist and can be modified if additional assert functions are needed.
  * There must be a central program called `run_tests.f90` which contains all the test calls defined in the modules.
  * Each subroutine's tests should be placed in independent modules (one file per tested subroutine). 
  * All test modules must be named `mod_<subroutine_name>.f90` to ensure they are compiled before `run_tests.f90`. Otherwise, compilation errors may occur.
  * Check details in `test/readme.md`

* **`/helper`** this folder will not be included in the final version of TOX. For now, it serves to help us create the C wrapper for the subroutines more quickly and easily. See details in `helper/readme.md`.

---

### FORD (Fortran Online Reference Documentation)

[FORD](https://github.com/Fortran-FOSS-Programmers/ford) is a documentation generator specifically designed for Fortran projects. It allows developers to create clean, structured, and navigable HTML documentation from source code using lightweight markup embedded in comments.

* Designed specifically for Fortran (unlike Doxygen which is general-purpose).
* Supports documentation of modules, subroutines, functions, derived types, and more.
* Uses `!!!` or `!>` comment syntax to annotate code.
* Ideal for scientific and engineering projects using modern Fortran.
* Easy to integrate into Git-based workflows.

Example usage:

```bash
ford ford.yml
```

This generates an HTML site you can explore in a browser (`doc/index.html` by default).


### Tensor Omics Snippets

Organize and place the snippets inside the appropriate snippet folders according to their functionality:

- Use the `f42:` prefix for F42-compliant infrastructure.
- Use the `tox:` prefix for application-specific Tensor Omics subroutines.

See `snippets/readme.md` for details.

---

# Compilation

> **Note:** If you're using prebuilt binaries from [Releases](https://gitlab.rlp.net/a.hallab/tensor-omics/-/releases), you can skip compilation and proceed directly to [Installation](#installation) → [Installation for R and Python](#installation-for-r-and-python) → [Usage](#usage).
>
> **If compiling from source:** Follow this workflow: Compilation → Installation → Installation for R and Python → Usage

Tensor Omics can be compiled using a **local toolchain** or via the reproducible **Docker environment**.

## Compiler Requirements

Tensor Omics requires **gfortran ≥ 15** for improved performance and modern Fortran support.

### Install Latest gfortran on Ubuntu/Debian

To ensure optimal performance and compatibility, we recommend adding a PPA (Personal Package Archive) to your system to easily install the latest version.

If you are using Ubuntu, you can do this by running:

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gfortran-15
```

Verify the installation:
```bash
gfortran-15 --version
```

### Alternative: Use Docker for Latest gfortran

If gfortran-15 is not available in your system repositories, you can also get the latest version with Docker:

- **Install and setup Docker** as explained for your operating system in the [Docker documentation](https://docs.docker.com/get-docker/).
- **Use our Dockerfile** `gfortran.docker` located in the `misc` directory:

```bash
cd misc
docker build -t arch-gfortran -f gfortran.docker .
cd ..
```

This creates a Docker image with the latest gfortran pre-installed. You can then use it for compilation (see [Compilation Using Docker](#compilation-using-docker) section for detailed instructions).



## Local Compilation (No Docker)

The `build.sh` script will compile all the files located in the `src/` directory.

It creates a directory for the compiled objects under `/build/<compiler>/`, and the resulting shared library will be named `libtensor-omics.so`.

This `.so` file is the one that must be loaded from R or Python.

Every time the code is compiled, a new `/build/<compiler>/` directory is created. To simplify access, the script creates a symbolic link to the latest compiled shared library so that R and Python can always load the same file consistently.

**Default compilation** (uses `gfortran` without performance optimizations):

```bash
./build.sh
```

**Maximum performance compilation** (uses `gfortran` with optimization flags):

```bash
./build.sh --max-performance
```

**Intel Fortran Compiler** (with maximum performance flags):

```bash
./build.sh --max-performance FC=ifx
```

> **Note:** Files are compiled in alphabetical order; please name your files accordingly.

## Compilation Using Docker

Install and setup Docker as explained for your operating system in the [Docker documentation](https://docs.docker.com/get-docker/).

### Step 1: Build the Docker Image

Navigate to the `misc` directory and build the Docker image:

```bash
cd misc
docker build -t arch-gfortran -f gfortran.docker .
```
### Step 2: Enter the Docker Container

Go back to the repository root and start an interactive Docker session:

```bash
cd ..
docker run -it -v $(pwd):/opt -w /opt arch-gfortran /bin/bash
```

### Step 3: Compile the Project

Inside the container, run:

```bash
./build.sh
```

### Step 4: Run Tests (Optional)

Inside the container, run:

```bash
./test_runner.sh
```

Or run tests for a specific module:

```bash
./test_runner.sh get_outliers
```


# Installation

Tensor Omics can be installed from prebuilt binaries or compiled from source.

## Prebuilt Binaries

Prebuilt binaries for Linux, Windows-11, and macOS are automatically generated by GitLab CI/CD pipelines and available in the [Releases](https://gitlab.rlp.net/a.hallab/tensor-omics/-/releases) section of this repository.

### Linux

1. Download the latest `.so` file from [Releases](https://gitlab.rlp.net/a.hallab/tensor-omics/-/releases)
2. Place it in a directory where R and Python can find it (e.g., `/usr/local/lib/` or your project directory)
3. Export the library path:
   ```bash
   export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
   ```

### Windows-11

1. Download the latest `.dll` file from [Releases](https://gitlab.rlp.net/a.hallab/tensor-omics/-/releases)
2. Place it in your system PATH or project directory
3. Ensure MSVC or GFortran runtime is installed:
   ```bash
   choco install mingw  # or vcredist for MSVC
   ```
4. Set the library path environment variable:
   ```cmd
   set PATH=C:\path\to\build;%PATH%
   ```
   Or for permanent setup, add to Windows Environment Variables.

### macOS

1. Download the latest `.dylib` file from [Releases](https://gitlab.rlp.net/a.hallab/tensor-omics/-/releases)
2. Place it in `/usr/local/lib/` or your project directory
3. Update the dynamic library search path:
   ```bash
   export DYLD_LIBRARY_PATH=/usr/local/lib:$DYLD_LIBRARY_PATH
   ```

---

# Installation for R and Python

Once the Tensor Omics shared library is available (either from prebuilt binaries or after compilation), follow these steps to use it from R and Python.

## R Integration

### Requirements

- **R ≥ 3.6**

### Installation Steps

1. **Load the library and source the wrapper functions**:
   ```r
   dyn.load("build/libtensor-omics.so")
   source("r/tensoromics_functions.R")
   source("r/error_handling.R")
   ```

2. **Verify installation**:
   ```r
   # Test a simple function
   result <- tox_euclidean_distance(c(1, 2, 3), c(4, 5, 6))
   print(result)  # Should output: 5.196152
   ```

### Example Usage

```r
# Load the library and wrapper functions
dyn.load("build/libtensor-omics.so")
source("r/tensoromics_functions.R")
source("r/error_handling.R")

# Use Tensor Omics functions
distance <- tox_euclidean_distance(c(0, 0), c(3, 4))
print(paste("Distance:", distance))
```

## Python Integration

### Requirements

- **Python ≥ 3.7**
- **numpy** (for array operations)

### Installation Steps

1. **Install required Python packages**:
   ```bash
   pip install numpy
   ```

2. **Load the library and import Tensor Omics functions**:
   ```python
   import numpy as np
   import ctypes
   import os
   from tensoromics_functions import euclidean_distance
   
   # Load library
   dll_path = os.path.abspath("build/libtensor-omics.so")
   ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
   lib = ctypes.CDLL(dll_path)
    ```

3. **Verify installation**:
    ```python
    distance = euclidean_distance([0.0, 0.0], [3.0, 4.0])
    print(f"Distance: {distance}")  # Should output: 5.0
    ```
### Example Usage

```python
import numpy as np
import ctypes
import os
from tensoromics_functions import (
    normalize_by_std_dev,
    quantile_normalization,
    euclidean_distance
)

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)

# Load gene expression data
expression_matrix = np.array([[1, 2], [3, 4], [5, 6]], dtype=np.float64)

# Apply normalization pipeline
normalized = normalize_by_std_dev(expression_matrix)
normalized = quantile_normalization(normalized)

# Calculate distances
distance = euclidean_distance(
    expression_matrix[:, 0],
    expression_matrix[:, 1]
)
print(f"Distance: {distance}")
```

---

# Usage



# Testing

The test suite framework provides a robust and scalable system for organizing and executing unit tests in Fortran.

## Test Architecture

1. **`run_tests.f90`** - Main program that handles command line arguments and runs all tests
2. **Test Modules** - Each module (suite) contains tests for a specific functionality
3. **`asserts.f90`** - Assertion function library for validating test results

## Running Tests Locally

```bash
# Run all tests from all suites
./test_runner.sh

# Run all tests from a specific suite
./test_runner.sh get_outliers

# Run specific tests from a suite
./test_runner.sh get_outliers test_compute_rdi_normal,test_compute_rdi_zero_scaling
```

## Running Tests in Docker

```bash
# Build and run all tests
docker run -it -v $(pwd):/opt -w /opt arch-gfortran ./test_runner.sh

# Run tests for a specific module
docker run -it -v $(pwd):/opt -w /opt arch-gfortran ./test_runner.sh get_outliers
```

## R and Python Integration Tests

After installation, verify R and Python integration:

### R Tests

```bash
# Run R-level tests
Rscript rcpp/test/euclidean_distance.R
Rscript rcpp/test/normalization.R
```

### Python Tests

```bash
# Run Python-level tests (if available)
python python/test/test_euclidean_distance.py
```


---

# Configuration

## Environment Variables

Set these variables to customize Tensor Omics behavior:

```bash
# Specify the number of OpenMP threads
export OMP_NUM_THREADS=8

# Set the library search path (Linux)
export LD_LIBRARY_PATH=/path/to/build:$LD_LIBRARY_PATH

# Set the library search path (macOS)
export DYLD_LIBRARY_PATH=/path/to/build:$DYLD_LIBRARY_PATH
```

**Windows:**
```cmd
# Add Tensor Omics shared library to PATH (temporary - current session only)
set PATH=C:\path\to\build;%PATH%

# Or for permanent setup, use setx (requires reopening terminal)
setx PATH "C:\path\to\build;%PATH%"
```

> **Note:** `/path/to/build` or `C:\path\to\build` must point to the directory containing:
> - `libtensor-omics.so` (Linux)
> - `tensoromics.dll` (Windows)
> - `libtensor-omics.dylib` (macOS)
## Build Options

Customize compilation via environment variables in `build.sh`:

| Variable | Default | Purpose |
|----------|---------|---------|
| `FC` | `gfortran` | Fortran compiler |
| `FFLAGS` | `-O3 -march=native` | Fortran compiler flags |


## Documentation Generation

Generate FORD documentation:

```bash
ford ford.yml
open doc/index.html
```

For details on documentation, see the [FORD documentation](https://github.com/Fortran-FOSS-Programmers/ford).

---

# License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.
