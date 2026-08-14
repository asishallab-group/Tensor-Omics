
# Tensor Omics

## About

Tensor Omics is a high-performance framework for explainable, geometry-based analysis of multimodal omics and related high-dimensional datasets. Instead of relying on black-box models, it treats expression profiles, clinical measures, or socioeconomic indicators as vectors in semantically meaningful spaces (e.g. tissues, disease stages, conditions). By measuring distances, angles, projections, and trajectories in these spaces, Tensor Omics enables direct comparison of activity across genes, paralogs, sexes, species, or patient groups. This geometric approach makes complex multivariate patterns interpretable and reproducible while remaining robust to sparsity and noise.

Designed for distributed high-performance computing, Tensor Omics is implemented in Fortran and C, using data-oriented parallel kernels and compiler-driven optimization (including vectorization and multicore execution when enabled). This makes the algorithms naturally decomposable into largely independent tasks, suitable for federated datasets where privacy and efficiency are critical. Scientific use cases include: detecting disease biomarkers and subtype-specific trajectories in medical data; quantifying divergence and neofunctionalization of gene duplicates in plant and animal transcriptomes; and reconstructing global gender-equality trajectories from socioeconomic indicators. Across these domains, Tensor Omics provides a unified, geometry-driven methodology for discovering explanatory patterns in heterogeneous, high-dimensional data.

### Key Features

* **Geometry-based analysis**: distances, angles, projections, and trajectory shifts are used as primary primitives.
* **Explainable outputs**: results are interpretable in terms of vector geometry rather than opaque model coefficients.
* **Multi-modal integration**: unifies transcriptomics, proteomics, metabolomics, clinical, or socioeconomic data within one framework.
* **Parallel and federated**: implemented in Fortran/C with data-oriented parallel kernels and compiler-level optimization for efficient large-scale computation on distributed datasets.
* **Robust to sparsity and noise**: percentile-based empirical thresholds and local geometric measures enable stability.
* **Broad applications**: demonstrated on medical biomarker discovery, gene duplication outcomes, developmental trajectories, and socioeconomic indicators.

*Tensor Omics shows that geometry, when treated not as preprocessing but as the central instrument of analysis, can open entirely new ways to read complex biological and social data — simple, transparent, and surprisingly powerful.*

Please read the manuscript `./misc/Tensor_Omics_Methods.pdf` for a detailed method description.

The **user manual** is currently available in the [`general_documentation`](https://github.com/asishallab-group/Tensor-Omics/tree/general_documentation) branch at `misc/tox_manual.pdf`. To access it locally:

```bash
git fetch origin general_documentation
git checkout general_documentation
# open misc/tox_manual.pdf with your preferred PDF viewer
```

### Project Structure

This repository contains the source code, methods, snippets, and tests for the **Tensor Omics (TOX)** project.

```
/build
  └── ...       # Compiled Fortran binaries and intermediate build files

/doc
  └── ...       # Documentation generated automatically using FORD (Fortran documentation tool)

/external
  └── ...       # Third-party Fortran sources (e.g. LOESS) and compiled external libraries

/material
  └── ...       # Example/reference datasets (.tsv) used by tests and demos

/misc
  └── ...       # Tensor Omics documentation, coding guides and helper dockerfile to compile the project

/python
  └── ...       # Python scripts that execute pipeline logic and invoke subroutines

/rcpp
  └── ...       # Rcpp scripts that execute pipeline logic and invoke subroutines

/snippets
  └── ...       # Code templates or reusable short logic blocks

/src
  └── ...       # Fortran backend

/test
  └── ...       # Fortran testing

/helper
  ├── helper_c_wrapper.py       # generates C wrapper subroutines
  ├── helper_rcpp_wrapper.py    # generates Rcpp wrapper subroutines
  └── generate_snippets.py      # auto-generates VS Code snippet files from source code

build.sh            # Compile and generate shared libraries (incl. bundled externals)
build_utils.sh      # Shared helper functions used by the build and test scripts
fpm.toml            # fpm compilation settings, also holds the FORD config
authors.h           # Author metadata macros used in the FORD documentation
test_runner.sh      # Compile and run Fortran unit tests
run_all_tests.sh    # Run the full test suite (Fortran + Python + R)
```

* **`/build`** stores shared libraries, compiled objects, and binary files resulting from Fortran compilation. It keeps the repository clean by separating source from compiled code.
* **`/doc`** contains auto-generated documentation built using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) from annotated Fortran source files.
* **`/external`** bundles third-party Fortran sources (e.g. the LOESS routines in `external/loess_netlib/`) that are compiled by `build.sh` into static libraries under `external/`.
* **`/material`** holds example and reference datasets (`.tsv`) used by the tests and usage demonstrations.
* **`/misc`** contains the team's coding guidelines (`misc/Fortran_Coding_Guides.pdf`), the detailed method description (`misc/Tensor_Omics_Methods.pdf`), and the Dockerfile (`misc/gfortran.docker`) used to compile the project without installing anything beyond Docker.
* **`/python`** includes Python scripts that coordinate analysis workflows.
* **`/rcpp`** includes Rcpp scripts that coordinate analysis workflows.
* **`/snippets`** includes frequently used or testable units of logic reused across development stages. Snippets expose subroutine names and their arguments. Use the `f42:` prefix for F42-compliant infrastructure and the `tox:` prefix for application-specific subroutines. See `snippets/readme.md` for details.
* **`/src`** contains performance-critical Fortran code compiled during the build process. Subroutines that perform no I/O operations or memory allocations must be declared `pure`.
* **`/test`** contains unit tests for the Fortran subroutines. `asserts.f90` provides the assertion library; `run_tests.f90` is the central test program. Test files (`mod_test_*.f90`) generally correspond to full modules from `src/`, though newly added subroutines or utility functions from utils are sometimes isolated into their own test cases. See `test/readme.md` for details.
* **`/helper`** is a temporary development aid for generating C wrappers (`helper_c_wrapper.py`), Rcpp wrappers (`helper_rcpp_wrapper.py`), and VS Code snippet files (`generate_snippets.py`). See `helper/readme.md` for details.

---

## Requirements

See [Install / Compilation](#install--compilation) for step-by-step setup instructions with Docker or native toolchain. For a native setup, you need:

| Requirement | Version | Notes |
|---|---|---|
| gfortran | ≥ 15 | Available natively on Arch Linux and macOS (Homebrew); use Docker on other systems |
| fpm | ≥ 0.12.0 | Fortran package manager |
| libzip | ≥ 1.11 | Required for archive handling |
| xxhash (XXH3) | ≥ 0.8 | Required for hashing |
| Python | ≥ 3.7 + NumPy + pandas | Python integration |
| R | ≥ 3.6 + Rcpp | R integration |
| gdb | optional | Needed only for `--debug` flag |

[FORD](https://github.com/Fortran-FOSS-Programmers/ford) (`pip install ford`) is optional and only needed to regenerate the API documentation.

---

## Get the code

Clone the repository from GitHub:

```bash
git clone https://github.com/asishallab-group/Tensor-Omics.git
cd Tensor-Omics
```

---

## Install / Compilation

We recommend using the provided **Docker image**: it comes with gfortran 15 and all required dependencies pre-installed, so no system-level changes are needed. If you already have gfortran 15+ installed, you can also compile natively.

### Using Docker (recommended)

Install Docker as explained for your operating system in the [Docker documentation](https://docs.docker.com/get-docker/).

**Step 1:** Build the Docker image from the Dockerfile provided in `misc/`:

```bash
cd misc
docker build -t arch-gfortran -f gfortran.docker .
cd ..
```

**Step 2:** Compile the project inside the container:

```bash
docker run -it -v $(pwd):/opt -w /opt arch-gfortran ./build.sh
```

### Native compilation

If you have **gfortran ≥ 15** installed, compile directly with `build.sh`. It compiles all files in `src/`, places compiled objects under `build/<compiler>_<hash>/`, and copies the resulting shared library (`libtensor-omics.so`) into `build/` so that Python and R always find it at the same path.

gfortran 15 is readily available on rolling-release or well-equipped systems. For most other Linux distributions and Windows, **we recommend Docker** since gfortran 15 is not yet available in standard package repositories (May 2026).

**Arch Linux** (rolling release — always has the latest GCC):

```bash
sudo pacman -S gcc-fortran
```

**macOS** (via [Homebrew](https://brew.sh)):

```bash
brew install gcc
```

On Ubuntu/Debian, gfortran 15 is not available in the standard repositories. Use the [Docker path](#using-docker-recommended) instead.

**Default** — gfortran with baseline optimisation flags:

```bash
./build.sh
```

**Maximum performance** — enables the `optimization` build profile (adds `-O3` and performance-oriented code paths):

```bash
./build.sh --max-performance
```

**Debug build** — enables debug/checking flags (`-O0`, debug symbols, and runtime checks depending on compiler profile):

```bash
./build.sh --debug
```

**Intel Fortran Compiler** — with additional optimisation flags:

```bash
FC=ifx ./build.sh --max-performance 
```

### Build options

Beyond the profiles above, `build.sh` accepts several options (which `test_runner.sh` inherits):

* `--compiler=<gfortran|ifx|nvfortran>` — compiler to use (defaults to `gfortran`). The `FC` environment variable is also honoured, with precedence `--compiler` > `$TOX_COMPILER` > `$FC`.
* `--max-performance` — enable the `optimization` profile (`-O3` and performance-oriented code paths).
* `--diagnostics` — enable diagnostic/debugging flags (helpful when debugging). Can be combined with `--max-performance`.
* `--override-flags="<flags>"` — replace the profile flags with your own, e.g. `--override-flags="-O2 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize -fPIC"`. When set, `--max-performance` has no effect.
* `--directive=<NAME>` — define a preprocessor directive; repeatable, e.g. `--directive=MAX_PERFORMANCE --directive=OTHER_DIRECTIVE`.
* `--clean-build` — force `fpm` to rebuild `src/` from scratch (enabled automatically when switching git branches). Useful when `fpm` misses changes that do not alter the module structure.

> **Note:** Each `--<option>` maps to an uppercased, `TOX_`-prefixed variable with non-alphanumeric characters replaced by underscores — e.g. `--override-flags` becomes `TOX_OVERRIDE_FLAGS`. Passing `--<option>=<value>` sets that value (a bare flag sets `1`), so any option can equivalently be supplied as an environment variable. An explicit `--<option>` always overrides the corresponding variable.

### Python integration

The `python/` folder contains the wrapper functions needed to call Tensor Omics subroutines from Python, along with example scripts in `python/test/` that demonstrate usage and can serve as a starting point for your own analyses.

Once the shared library is available, load it from Python:

```python
import ctypes, os

lib_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)

lib = ctypes.CDLL(lib_path)
from tensoromics_functions import tox_euclidean_distance
```

Verify that the library loaded correctly:

```bash
python3 python/test/mod_test_euclidean_distance.py
```

### R integration

The `rcpp/` folder contains the Rcpp wrapper functions for calling Tensor Omics from R, along with example tests in `rcpp/test/` that illustrate how to use each function.

```r
library(Rcpp)
lib_path <- shQuote(normalizePath("build"))
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))
Rcpp::sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)
```

Verify that the library loaded correctly:

```bash
Rscript rcpp/test/mod_test_euclidean_distance.R
```

---

## API Documentation and Snippets

### Reading the Fortran API documentation

The `doc/` folder contains pre-generated HTML documentation produced by [FORD](https://github.com/Fortran-FOSS-Programmers/ford) from the annotated Fortran source files. Open it in any browser:

```bash
# Linux
xdg-open doc/index.html

# macOS
open doc/index.html
```

The documentation covers all public modules, subroutines, functions, and derived types in `src/`. Navigate through the sidebar to browse by module, or use the search bar to look up a specific subroutine by name.

To regenerate the documentation after making changes to the source code, install FORD (`pip install ford`) and run:

```bash
ford concise_project_info.md
```

### Loading snippets into VS Code

The `snippets/` folder contains VS Code snippet files for Fortran, Python, and R. Loading them gives you prefix-triggered autocomplete for all Tensor Omics subroutines — type a prefix and press `Tab` to insert the full call with argument placeholders.

Snippets are organized by prefix:

* `tox:*` — Tensor Omics analysis subroutines (start here if you are doing analysis)
* `f42:*` — general-purpose infrastructure subroutines (sorting, data structures, etc.)
* `toxdev:*` — development convenience snippets (parallel loop templates, argument declarations, etc.)

**To load snippets in VS Code:**

1. Open the Command Palette (`Ctrl+Shift+P` / `Cmd+Shift+P`).
2. Type `Snippets: Configure Snippets` and press Enter.
3. Select the target language (e.g. `fortran`, `python`, `r`) or choose `New Global Snippets file...` to make all snippets available in any file.
4. Copy the contents of the corresponding file from `snippets/` (e.g. `Fortran_tox_snippets.json` for Fortran analysis snippets) and paste it into the snippet file that opened.

Alternatively, copy the snippet files directly to your VS Code user snippets directory:

```bash
# Linux
cp snippets/Fortran_tox_snippets.json ~/.config/Code/User/snippets/
cp snippets/Python_tox_snippets.json ~/.config/Code/User/snippets/
cp snippets/R_tox_snippets.json       ~/.config/Code/User/snippets/

# macOS
cp snippets/Fortran_tox_snippets.json "$HOME/Library/Application Support/Code/User/snippets/"
cp snippets/Python_tox_snippets.json  "$HOME/Library/Application Support/Code/User/snippets/"
cp snippets/R_tox_snippets.json       "$HOME/Library/Application Support/Code/User/snippets/"
```

Restart VS Code after copying. Snippets are auto-generated from the source code by the CI pipeline — see `snippets/readme.md` for details on how they are maintained.

---

## Testing

The test suite provides a scalable system for organizing and executing Fortran unit tests.

### Architecture

1. **`run_tests.f90`** — main program handling command-line arguments and dispatching all test calls.
2. **`test_suite.f90`** — The core orchestration module. It defines abstract interfaces, registries, and runner subroutines to dynamically collect and execute tests across different suites.
3. **Test modules** — Test files that generally correspond to complete modules in `src/`. They may also be dedicated to specific standalone or newly added subroutines, particularly those from utils.
4. **`asserts.f90`** — assertion library for validating results.

* **Suite Registration:** Each test module must implement a public function named `get_all_<src_module_name>_tests` that matches the `get_all_interface` abstract interface defined in `test_suite.f90`. This function aggregates and returns an array of all `test_case` structures defined within that module.
* Test modules must be named `mod_test_<module_name>.[fF]90` to ensure they are compiled before `run_tests.f90`. Files are compiled in alphabetical order; name test files accordingly. See `test/readme.md` for details.

### Running tests

```bash
# Run the full test suite (Fortran + Python + R)
./run_all_tests.sh

# Run only Fortran tests
./test_runner.sh

# Run all tests from a specific Fortran suite
./test_runner.sh <suite_name>

# Run specific tests within a suite
./test_runner.sh <suite_name> <test1,test2,test3>
```

### Debugging tests

If tests fail or you need to inspect variables at runtime, run the Fortran test runner in debug mode:

```bash
./test_runner.sh --debug
```

This launches the test executable under `gdb`, so you can set breakpoints, inspect variables and call stacks, and step through execution.

You can also run specific suites or tests in debug mode:

```bash
./test_runner.sh --debug <suite_name>
./test_runner.sh --debug <suite_name> <test1,test2,test3>
```

If `gdb` is not installed on your system, install it first. For command reference, see the [gdb manual](https://sourceware.org/gdb/current/onlinedocs/gdb/).

### Test runner options

`test_runner.sh` accepts all of the [build options](#build-options) above, plus a few test-specific ones:

* `--skip-kinds-test` — skip the compile-time check that C types match the Fortran kinds. Handy to avoid the extra clean build it otherwise triggers, since the check always passes on an unchanged platform.
* `--reuse-mod-files` — keep `fpm`'s test module files instead of removing them before each run. Speeds up test recompilation; use when debugging tests.
* `--test-target=<target>` — select the test target from `fpm.toml` (currently only `run_tests`, which is the default).
* `--keep-files` — keep the temporary files the runner creates in the repo root (removed by default).
* `--keep-<ext>` — the fine-grained variant of `--keep-files`, keeping only files of a given extension, e.g. `--keep-zip` or `--keep-txt`.

---

## Contribute

Tensor Omics is developed openly and contributions are welcome — whether reporting issues, improving documentation, or suggesting new analysis subroutines.

When writing Fortran code, consult the team's coding guidelines in `misc/Fortran_Coding_Guides.pdf`. 
API documentation is generated from annotated source comments — see [Reading the Fortran API documentation](#reading-the-fortran-api-documentation) for how to regenerate it.

