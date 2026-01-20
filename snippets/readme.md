### How to Add Snippets in VSCode

To install a snippet in Visual Studio Code:

1. Open the command palette (`Ctrl+Shift+P` or `Cmd+Shift+P`).
2. Search for:  
   `Preferences: Configure Snippets`
3. Choose a global file (e.g. `global.code-snippets`) or the language-specific file for Fortran (e.g. `fortran.json`).
4. Paste the snippet block inside the JSON file.

---

### Snippets explanation

There are three kinds of snippets:
1. `f42:*`: These snippets are related to functions/subroutines that serve general purposes, like sorting an array or creating a kd tree.
2. `tox:*`: These snippets are related to anything TensorOmics-related, so things that are needed for analysises
3. `toxdev:*`: These snippets are just convenience snippets for development, like the below `do-parallel` snippet, or snippets for argument declaration.

#### `toxdev:do-parallel`

Start typing `toxdev:do-parallel` and let the snippet autocomplete.

It includes tab stops for easy customization:
- `${1:i}` → loop index variable
- `${2:N}` → name of the loop limit variable
- `${3:10}` → initial value assigned to `N`
- Loop body → the cursor jumps there after filling in the loop setup

The snippet expands into a unified parallel loop structure that automatically selects the appropriate backend based on compiler flags:
- Coarray
- OpenMP
- Serial fallback

This allows backend-independent parallel execution using preprocessor directives.

---

### Developer Notes

The Snippets are being generated automatically from the repository's code by the CI pipeline.

The corresponding script for generation is located at `helper/generate_snippets.py`

To make this work, there are some requirements (otherwise helpful errors are thrown):
1. Each declared function/code in Python/R (not Fortran) needs a preceding header comment, with patterns:
    - For functions related to real Fortran functions: `#> <fortran_module_name>:<fortran_wrapper_name>:<description>`: 
      <details><summary>See Examples</summary>

      ```python
      #> f42_serialize_real:serialize_real_nd_C: Serialize an n-dimensional array of type 'float'
      def tox_serialize_real_nd(arr: np.ndarray, filename: str):
          """
          Serializes an n-dimensional real64 array to a binary file
          """
          ...
      ```
      ```r
      #> tox_euclidean_distance:euclidean_distance_c: Calculate Euclidean distance between two vectors
      #' 
      #' Computes the Euclidean distance between two vectors of the same dimension.
      #' This function automatically checks for errors and throws informative exceptions.
      #' 
      #' ...
      tox_euclidean_distance <- function(vec1, vec2) {
        ...
      }
      ```
      </details>

    - For helper functions that don't wrap a specific Fortran routine: `#> (f42|tox)_helper:<description>`:
      <details><summary>See Examples</summary>

      ```python
      #> f42_helper: convert a filename to ASCII chars to transfer it as integer to fortran
      def _filename_to_ascii_array(filename):
       ascii_arr = np.array([ord(c) for c in filename], dtype=np.int32)
       return ascii_arr, np.int32(len(ascii_arr))
      ```
      ```python
      #> f42_helper: Alias for build_kd_index to build a spherical KD-Tree index for the given unit vectors
      def build_spherical_kd(vectors, dimension_order=None):
       """
       Alias for build_kd_index to build a spherical KD-Tree index for the given unit vectors.
      
       Parameters:
       vectors (np.array): 2D array of unit vectors (d x n, Fortran order)
       dimension_order (np.array): Order of dimensions for splitting (1-based)
      
       Returns:
       np.array: Spherical KD-Tree indices (1-based Fortran indices)
       """
       # For spherical KD-Tree, we use the same implementation as regular KD-Tree
       # but with a different name for clarity
       return build_kd_index(vectors, dimension_order)
      ```
      ```r
      #> f42_helper: Clean data by removing or imputing problematic values
      #' Clean data by removing or imputing problematic values
      #'
      #' This function handles NA, NaN, Inf values and genes that are all zeros
      #' to prepare data for Fortran normalization routines.
      #' ...
      tox_clean_data_for_normalization <- function(df_matrix, 
                                           remove_all_zero_genes = TRUE,
                                           na_strategy = "remove_genes",
                                           min_expression_threshold = 0.0,
                                           convert_small_to_zero = FALSE) {
      ...
      }
      ```
      </details>

    - For code without a declared function: `#> (f42|tox)_helper-<name>:<description>`:
      <details><summary>See Examples</summary>

      ```python
      #> f42_helper-import_libs: Import necessary packages
      import numpy as np
      import ctypes
      import os
      ```
      ```r
      #> f42_helper-import_libs: Import necessary packages
      library(Rcpp)
      
      # Get absolute path to build directory containing the compiled Fortran library
      
      lib_path <- shQuote(normalizePath("build"))
      
      # Set up compilation flags for linking with Fortran library
      Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))
      
      # Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
      sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)
      
      cat("✓ TensorOmics Rcpp functions loaded successfully\n")
      
      source("rcpp/error_handling.R")
      ```
      </details>

    - **_Note_** that all empty lines or comments after the last code line of a snippet will be removed:
      <details><summary>See Examples</summary>

      ```python
      #> f42_helper-import_libs: Import necessary packages
      import numpy as np
      import ctypes
      import os

      # everything from here will be removed, as it is either comment
      # bla
      """ or empty line """

      # to avoid garbage and allow section comments like
      # -----------------------------------
      # Fancy Stuff from Here
      # -----------------------------------
      ```
      ```r
      #> f42_helper-import_libs: Import necessary packages
      library(Rcpp)
      
      # Get absolute path to build directory containing the compiled Fortran library
      
      lib_path <- shQuote(normalizePath("build"))
      
      # Set up compilation flags for linking with Fortran library
      Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))
      
      # Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
      sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)
      
      cat("✓ TensorOmics Rcpp functions loaded successfully\n")
      
      source("rcpp/error_handling.R")
      
      # From here everything is ignored
      # ===================================================================
      # EUCLIDEAN DISTANCE FUNCTIONS
      # ===================================================================
      ```
      </details>
    - **__Note__** that all lines before the first appearing line starting with `#>` in a Python or R file will be ignored.

2. The snippet generation expects valid Python/R/Fortran syntax, this is not a linter.
3. This generator enforces the anyway proposed module naming pattern `(f42|tox)_.*`, so each module name needs a matching prefix for what it is related to, either `tox` or `f42`.
4. Each defined Fortran wrapper needs a defined Python function. E.g. if there is `subroutine_c` for module `f42_module` in Fortran, the generation will fail if there is no related `#> f42_module:subroutine_c: Bla bla` in Python. For R this is not checked because of the current switch from `.Fortran` to `Rcpp`
5. There is a special header implemented to skip all lines below it. The line must start with `#>skip snippets` in this case. Anyway, this shouldn't be needed. It was only added to allow the error handling in R to have the `check_err_code` at the beginning, while defining input validation functions afterwards that don't need snippets.
6. The following files will be included for snippet generation:
    - `python/*.py`
    - `rcpp/*.[rR]`
    - `src/**/*.[fF]90`, except `src/config.F90` and `src/safeguard.F90`
    - In case that more paths need to be included/excluded, just adjust the `main` function of `generate_snippets.py`. All `generate_*_snippets` functions support an `ignored_files: List[str]` argument, a list of paths.
