# Subroutine Interface Generator

This Python script helps generate a **Fortran subroutine wrapper** that can be called from **C or Python** via shared libraries. It is useful within the **Tensor Omics** project to ensure consistent and reusable interfaces for  routines written in Fortran.

---

## What does this script do?

The script automatically generates a Fortran-compatible interface (`bind(C)`) for an existing subroutine, making it easier to:

- Interoperate with C or Python (via `ctypes`)
- Maintain a consistent call signature
- Reuse argument declarations cleanly

---

## What should you modify for each subroutine?

To adapt the script for a new core subroutine, you need to update the following variables:

### 1. `core_name`
The name of the main subroutine you want to expose.

```python
core_name = "quantile_normalization"
```

---

### 2. `module_name`

The name of the Fortran module where the subroutine is defined.

```python
module_name = "tox_normalization_mod"
```

---

### 3. `core_args`

A multiline string with the list of arguments exactly as declared in the Fortran source. 

```python
core_args = """
  integer, intent(in) :: n_genes, n_tissues, max_stack
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  ...
"""
```

---

### 4. `arg_order`

A list with the **exact order** in which arguments should be passed to the subroutine. This may differ from the order in `core_args`, especially if grouped by intent.

```python
arg_order = [
    "n_genes", "n_tissues", "input_matrix", "output_matrix",
    "temp_col", "rank_means", "perm", "stack_left", "stack_right", "max_stack"
]
```

---

## Output

The script generates a `bind(C)` Fortran wrapper with a standard C-compatible interface, which can be compiled into a shared object (`.so`) and loaded from C or Python.

It returns a string in the terminal that you can copy and paste into your wrapper file.

