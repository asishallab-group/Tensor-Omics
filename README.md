
# Tensor Omics

Tensor Omics is a high-performance framework for explainable, geometry-based analysis of multimodal omics and related high-dimensional datasets. Instead of relying on black-box models, it treats expression profiles, clinical measures, or socioeconomic indicators as vectors in semantically meaningful spaces (e.g. tissues, disease stages, conditions). By measuring distances, angles, projections, and trajectories in these spaces, Tensor Omics enables direct comparison of activity across genes, paralogs, sexes, species, or patient groups. This geometric approach makes complex multivariate patterns interpretable and reproducible while remaining robust to sparsity and noise.

Designed for distributed high-performance computing, Tensor Omics is implemented in Fortran and C with OpenMP parallelisation, SIMD optimisation, and Fortran Coarrays, making the algorithms embarrassingly parallel and suitable for federated datasets where privacy and efficiency are critical. Scientific use cases include: detecting disease biomarkers and subtype-specific trajectories in medical data; quantifying divergence and neofunctionalization of gene duplicates in plant and animal transcriptomes; and reconstructing global gender-equality trajectories from socioeconomic indicators. Across these domains, Tensor Omics provides a unified, geometry-driven methodology for discovering explanatory patterns in heterogeneous, high-dimensional data.

## Detailed Method Description

Please read the manuscript `./misc/Tensor_Omics_Methods.pdf` for details.

## Key Features

* **Geometry-based analysis**: distances, angles, projections, and trajectory shifts are used as primary primitives.
* **Explainable outputs**: results are interpretable in terms of vector geometry rather than opaque model coefficients.
* **Multi-modal integration**: unifies transcriptomics, proteomics, metabolomics, clinical, or socioeconomic data within one framework.
* **Parallel and federated**: implemented in Fortran/C with OpenMP, SIMD, and coarrays for efficient large-scale computation on distributed datasets.
* **Robust to sparsity and noise**: percentile-based empirical thresholds and local geometric measures enable stability.
* **Broad applications**: demonstrated on medical biomarker discovery, gene duplication outcomes, developmental trajectories, and socioeconomic indicators.


*Tensor Omics shows that geometry, when treated not as preprocessing but as the central instrument of analysis, can open entirely new ways to read complex biological and social data — simple, transparent, and surprisingly powerful.*

---

# TOX Project Structure

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
└── ...         #Code templates or reusable short logic blocks 

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


### Compilation

The `build.sh` script will compile all the files located in the `src/` directory.

It creates a directory for the compiled objects under `/build/<compiler>/`, and the resulting shared library will be named `libtensor-omics.so`.

This `.so` file is the one that must be loaded from R or Python.

Every time the code is compiled, a new `/build/<compiler>/` directory is created. To simplify access, the script creates a symbolic link to the latest compiled shared library so that R and Python can always load the same file consistently.

Usage:

→ Uses the `gfortran` compiler without performance optimizations.

```bash
./build.sh
```

→ Uses the `gfortran` compiler with maximum performance flags.

```bash
./build.sh --max-performance
```

→ Uses the `ifx` compiler with maximum performance flags.

```bash
./build.sh --max-performance FC=ifx
```

Keep in mind that files are compiled in alphabetical order, please name your files accordingly.

---


### Testing

The test suite framework provides a robust and scalable system for organizing and executing unit tests in Fortran. It allows running individual tests, complete test suites, or all project tests with simple and clear syntax.

#### Architecture

1. **`run_tests.f90`** - Main program that handles command line arguments
2. **Test Modules** - Each module (suite) contains tests for a specific functionality
3. **`asserts.f90`** - Assertion function library for validating results

#### System Usage

```bash
# Run all tests from all suites
./test_runner.sh

# Run all tests from a specific suite
./test_runner.sh <suite_name>

# Run specific tests from a suite
./test_runner.sh <suite_name> <test1,test2,test3>
```

Keep in mind that files are compiled in alphabetical order, please name your files accordingly.

See `test/readme.md` for details.

---

## Get latest gfortran with Docker

Install and setup Docker as explained for your operating system in the Docker
documentation.

Use our Dockerfile `gfortran.docker` in `misc` directory:
```bash
docker build -t arch-gfortran -f gfortran.docker .
```

Then build the project with:
```bash
docker run -it -v `pwd`:/opt arch-gfortran ./build.sh
```

Use `./test_runner.sh` if you want to run the unit tests for the modules. In case you want to test only one module, use `./test_runner.sh <test suite name>`, e.g. `./test_runner.sh get_outliers`

Feel free to extend this README with additional information.
