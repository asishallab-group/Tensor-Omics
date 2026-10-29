Tensor Omics is a high-performance framework for explainable, geometry-based analysis of multimodal omics and related high-dimensional datasets. Instead of relying on black-box models, it treats expression profiles, clinical measures, or socioeconomic indicators as vectors in semantically meaningful spaces (e.g. tissues, disease stages, conditions). By measuring distances, angles, projections, and trajectories in these spaces, Tensor Omics enables direct comparison of activity across genes, paralogs, sexes, species, or patient groups. This geometric approach makes complex multivariate patterns interpretable and reproducible while remaining robust to sparsity and noise.

Designed for distributed high-performance computing, Tensor Omics is implemented in Fortran and C with OpenMP parallelisation, SIMD optimisation, and Fortran Coarrays, making the algorithms embarrassingly parallel and suitable for federated datasets where privacy and efficiency are critical. Scientific use cases include: detecting disease biomarkers and subtype-specific trajectories in medical data; quantifying divergence and neofunctionalization of gene duplicates in plant and animal transcriptomes; and reconstructing global gender-equality trajectories from socioeconomic indicators. Across these domains, Tensor Omics provides a unified, geometry-driven methodology for discovering explanatory patterns in heterogeneous, high-dimensional data.

## Key Features

* **Geometry-based analysis**: distances, angles, projections, and trajectory shifts are used as primary primitives.
* **Explainable outputs**: results are interpretable in terms of vector geometry rather than opaque model coefficients.
* **Multi-modal integration**: unifies transcriptomics, proteomics, metabolomics, clinical, or socioeconomic data within one framework.
* **Parallel and federated**: implemented in Fortran/C with OpenMP, SIMD, and coarrays for efficient large-scale computation on distributed datasets.
* **Robust to sparsity and noise**: percentile-based empirical thresholds and local geometric measures enable stability.
* **Broad applications**: demonstrated on medical biomarker discovery, gene duplication outcomes, developmental trajectories, and socioeconomic indicators.

## How the API is shaped

Every procedure has one name in Fortran, Python and R. From C you call it by that name with the suffix `_c` (`compute_rdi_c`): a C-compatible Fortran wrapper, kept in a module whose name ends in `_c` as well. That is why these pages list a module up to three times: `tox_get_outliers_impl` holds the implementation, `tox_get_outliers` is the Fortran API and `tox_get_outliers_c` the C API. A procedure whose name ends in `_impl` is never published, because it skips all input validation. A procedure's documentation, found on these pages, is its documentation in every language.

**Two tiers, where there is something to take over.** `foo` is the entry point to reach for first: it validates its arguments, allocates the work arrays and prepares what it can derive itself (a sorted permutation, a workspace size, a threshold), then computes. `foo_expert`, where it exists, validates and computes with what you hand it. It allocates nothing and prepares nothing, so you can reuse your own buffers or supply a different sort order. Fortran and C always get both tiers. Python and R get the expert tier only where it offers more than buffers: their `foo` calls the library's `foo`, which allocates its own work arrays, so an expert tier that only takes buffers would save them nothing.

**The calling language allocates the outputs.** No procedure allocates an array and hands it back. In Fortran and C you pass the output arrays in. In Python and R the binding allocates them as NumPy arrays or R vectors before the call, and the library writes into that memory directly. Every result is therefore an ordinary object of your language: its garbage collector frees it, and there is nothing to release by hand. The only memory the library allocates itself is scratch space inside `foo`, and it is released before `foo` returns.

This is deliberate. An all-in-one entry point that allocates its outputs inside the library and returns them would save the bindings a few lines. But memory the library allocated must also be freed by the library, so every language would need its own finaliser calling back into it, or a copy of each result into a native array. The cleanup would still be needed, only harder to get right and in more places.

### In every language

**Errors come back in `ierr`,** the last argument of every procedure. `0` is success. Anything else is `10000*arg_pos + code`: the last four digits are the code, which names the failure (1xx input and output, 2xx invalid input, 3xx memory, 5xxx Fortran runtime, 9xxx internal), and the digits before them are `arg_pos`, the argument to blame, counted in the Fortran signature and `0` when no argument is at fault. Python turns it into an exception derived from `ToxError` (`ToxInputError`, `ToxIOError`, ...) and R into a condition of class `tox_error` (`tox_input_error`, ...), both naming the argument. An argument Python or R refuses before the call, for its shape, type or value, is an input error too but carries no code: Python's `ToxArgumentError` (a `ToxInputError`, and also a `TypeError` and a `ValueError`) with `code` None, R's `tox_type_error`, `tox_shape_error` or `tox_na_error` (each a `tox_input_error`) with `code` NA. So `except ToxInputError` and `tryCatch(tox_input_error = ...)` catch every input error, the library's and the binding's alike. The codes are named: `ERR_*` for errors, `STAT_*` for outcomes that are not failures, such as an entry of a status array. Python exports them, every exception class, `ARG_POS_FACTOR` and `check_err_code` from the package itself (`from tensor_omics import ERR_DIVISION_BY_ZERO, ToxInputError`), and R loads them as globals.

**Indices are 1-based,** as in Fortran: family indices and permutations alike, in what you pass as much as in what comes back. Where an index can be "none", its value is `0`, such as a gene that belongs to no family.

**Trailing blanks are not part of a string,** as in Fortran: `"genes.tsv  "` and `"genes.tsv"` are the same value. Python and R drop them from what you pass, mode arguments included (`"ward "` selects `ward`), and a string never comes back to them with trailing blanks; from C, trim them yourself (below). A leading blank is part of the value.

**Results never contain NaN.** Where an output value is undefined, the documentation of that output names a sentinel outside the range the quantity can take, such as `-1` for a distance that was not measured. A NaN in a result is a bug; please report it. On input, NaN and infinite values are rejected unless an argument's documentation says they are permitted.

### From C

Link against `libtensor_omics.so`. No header is shipped; the prototypes follow from the `_c` modules on these pages.

- Every argument is passed by pointer. The types are `double`, `int` (32 bit), `bool` (`_Bool`) and `char`.
- Arrays are column-major, as in Fortran. Their extents are arguments: either named ones such as `n_genes`, or `n_<name>_elements` right after the array. Arrays that share an extent must agree on it; the C layer does not check that, Python and R do. Where the documentation declares an axis by an expression, such as `expr(sum(reps_per_tissue), n_genes)`, you pass that axis's actual extent as `n_expr_elements_dim_1` right after the array, and the C wrapper returns `ERR_SIZE_MISMATCH` unless it equals the expression.
- A string is a `char` buffer followed by its width, `<name>_strlen`. An array of strings is one column-major block of equally wide strings, followed by `<name>_strlen` and `n_<name>_elements`. **Pad with blanks, do not terminate with NUL**: `"genes.tsv"` in a `char[16]` passes a filename that ends in NULs. Fill an output buffer with blanks before the call and trim the trailing blanks from what comes back, so a value that itself ends in blanks does not survive the trip. Mode arguments (`mode`, `method`) are the exception and accept either convention.
- An optional argument without a default is absent when you pass `NULL`. An argument with a default is required in C, and its documentation states the default.
- `foo_expert_c` takes the workspace sizes that `foo_c` computes for you. The documentation of each size names the routine that computes it.
- An output whose valid length is known only after the call has its full size, and a count argument says how much of it is valid, such as `n_unique`. Python and R return only the valid part.

### From Python and R

Python finds the library through `TENSOR_OMICS_LIBRARY` or at `build/libtensor_omics.so`, and R loads it with `r/load_tensor_omics.R`. Neither asks for an extent that follows from the arrays you pass, nor for workspace sizes or string widths. An extent that only sizes an output, such as `n_families`, is still an argument. Python passes arrays to the library in column-major order, so a row-major NumPy array is copied first; pass `np.asfortranarray` data to avoid that for large inputs. Python results are read-only: `.copy()` one to modify it. An argument the procedure changes in place is changed in place in Python and returned as a new value in R.

**What an argument accepts** is the same in Python and R, and follows from its Fortran type, which each procedure's documentation shows. A value is accepted only if it converts to that type and back unchanged. Nothing is rounded, truncated, wrapped around or parsed on the way, and a refusal names the argument.

| Fortran type | Accepted | Refused |
|---|---|---|
| `integer(int32)` | Python and NumPy integers, and floats whose values are whole, as long as every value fits in 32 bits; R integers, and doubles with whole values (`2` is a double in R) | booleans and R logicals, fractions (`1.7`), NaN and NA, infinities, values out of range (`2**32 + 5`), strings, datetimes, complex numbers, R factors, Dates and `integer64` |
| `real(real64)` | any Python or NumPy integer or float that a float64 holds exactly, which is every float up to float64 and every integer up to `2**53` in magnitude; R integers and doubles. NaN and infinite values only where the argument's documentation permits them | booleans and R logicals, an integer or a long double a float64 would round (`2**53 + 1`; convert it yourself where the rounding is acceptable), strings (never parsed as numbers), complex numbers, datetimes, R Dates |
| `complex(real64)` | NumPy arrays of numbers, complex ones included, converted as NumPy converts them; R complex and numeric vectors. A single complex result comes back as a Python `complex` or an R complex | in R, anything that is not a number |
| `logical` | booleans, R logicals, and numbers that are all exactly `0` or `1`, such as `c(1, 0, 1)` | any other number, NaN, NA, `None`, strings |
| a string | a Python `str` or an `os.PathLike` path; an R character vector of length 1 | bytes, `None`, lists, numbers, NA, a NUL character, invalid UTF-8, and a value longer than a fixed declared length (never cut short) |
| an array of strings | a list, tuple, NumPy array or pandas Series of `str`; an R character vector, or a factor (its labels) | a bare `str`, a set or dict (they have no order), any element that is not a string, NA, a NUL character, invalid UTF-8, an element longer than a fixed declared length |
| a mode | one of the strings the documentation lists, matched whole, case ignored | anything else; the message lists the modes |

Arrays must have the declared rank: nothing is flattened, so a matrix passed for a vector is refused. In R a plain vector, a 1-D array and a table count as vectors, and an all-numeric data.frame as a matrix. A NumPy masked array is accepted only with no masked elements. R strings marked latin1 reach the library translated to UTF-8, and strings come back to R marked UTF-8. An R string with no encoding mark, such as one `readLines` gives, is passed as it stands if its bytes are UTF-8, in any locale; otherwise, and for a string marked "bytes", it is refused: mark a latin1 one with `Encoding(x) <- "latin1"`.

Where a refusal has an obvious fix, the message names it: `.view(np.int32)` for a `uint32` bit-mask word of 2^31 or more, `.astype(str)` for a bytes array, `.filled()` for a masked array, `as.integer()` or `as.double()` for an R `integer64`.


*Tensor Omics shows that geometry, when treated not as preprocessing but as the central instrument of analysis, can open entirely new ways to read complex biological and social data — simple, transparent, and surprisingly powerful.*
