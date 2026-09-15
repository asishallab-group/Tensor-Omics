# Tests

TensorOmics is tested in three languages, each with its own job (see "Where to Test What" in
`misc/Fortran_Coding_Guides.tex`):

- **Fortran** is the single source of truth for numerical correctness: expected values, every
  early exit, every bound.
- **Python and R** test the bindings only: every published procedure can be called, returns the
  documented type and shape, and raises the documented error. The bindings are generated, but their
  tests are not.

Every module has one suite per language, named `mod_test_<module>`.

## Layout

```
test_framework/                          # The Fortran test framework, a local fpm package
├── fpm.toml                             # (a dev-dependency of tensor_omics, see below)
└── src/
    ├── asserts.F90                      # Assertions
    └── test_suite.F90                   # Suite registry, running, reporting
test/
├── run_tests.F90                        # The test program: registers every suite
├── mod_test_sorting.F90                 # A suite in a single file
├── mod_test_tox_normalization/          # A suite with one child module per procedure
│   ├── mod_test_tox_normalization.F90   # The suite itself, gathering the children's cases
│   ├── mod_test_calc_fchange.F90
│   ├── ...
│   └── mod_test_tox_normalization_fixtures.F90   # Data the children share
└── test_files/                          # Input files for the tox_data suites
python/test/mod_test_<module>.py         # Python binding suites
r/test/mod_test_<module>.R               # R binding suites
```

The framework is a package of its own because fpm lets a test source use only library modules and
modules in its own directory or below. As a library, `asserts` and `test_suite` are visible to a
suite in any subdirectory of `test/`. It uses only the intrinsic modules: tensor_omics depends on
it, so it cannot depend back on tensor_omics.

## Running

```bash
./test_runner.sh                             # every Fortran suite
./test_runner.sh tox_normalization           # one suite
./test_runner.sh tox_normalization test_calc_fchange_values,test_log2_values   # some cases
./run_all_tests.sh                           # Fortran, then every Python and R suite
```

`test_runner.sh` builds the library and the tests, then runs `run_tests`. Before that, it runs the
kinds test: eight builds, each forcing one C-kind mismatch that `f42_safeguard` must reject. Set
`TOX_SKIP_KINDS_TEST=1` to skip it while iterating. `run_all_tests.sh` exits non-zero if any part
fails.

Both take the build options of `build.sh`. Every `--name[=value]` can also be set as the environment
variable `TOX_NAME`:

| Option | Effect |
|---|---|
| `--compiler=gfortran\|ifx\|nvfortran` | the compiler (or `$FC`, `$TOX_COMPILER`); gfortran by default |
| `--diagnostics` | runtime checks and every warning the compiler has (`-fcheck=all -Wall -Wextra ...`) |
| `--max-performance` | optimized build |
| `--debug` | run the tests under `gdb` |
| `--directive=NAME` | define the preprocessor symbol `NAME`, e.g. `NO_R_BINDING` |
| `--override-flags="..."` / `--override-link-flags="..."` | replace the compiler / link flags |
| `--clean-build` | rebuild everything |
| `--skip-code-generation` | build the generated bindings as they are in the tree |

The test runner has a few variables of its own: `TOX_TEST_TARGET` names another test program, and
`TOX_KEEP_FILES`, `TOX_KEEP_ZIP`, `TOX_KEEP_TXT` and `TOX_KEEP_BIN` keep the files the `tox_data`
suites write.

### Output

Each case prints `✓ <name> passed.` or `✗ <name> failed (<n> assertion(s)).`, followed by the
failed assertions. A failed assertion does not stop its case: every check in it runs and is
reported. The run ends with `All <n> test cases passed.`, or with how many failed across how many
suites, and then exits with status 1.

## Writing Fortran tests

### A suite in one file

```fortran
!> The `tox_foo` cases: ...
module mod_test_tox_foo
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_foo
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_tox_foo() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(2))
        all_tests(1) = test_case("test_foo_values", test_foo_values)
        all_tests(2) = test_case("test_foo_dimensions", test_foo_dimensions)
    end function get_all_tests_tox_foo

    !> What this case establishes, and how the expected values are derived.
    subroutine test_foo_values()
        real(real64) :: expr(2, 2), result(2, 2), expected(2, 2)
        integer(int32) :: ierr

        expr(:, 1) = [1.0_real64, 7.0_real64]
        ...
        call foo(2, 2, expr, result, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_foo_values: ierr")
        call assert_equal_array_real(result, expected, 4, 1e-12_real64, "test_foo_values: ...")
    end subroutine test_foo_values
    ...
end module mod_test_tox_foo
```

Register it in `test/run_tests.F90`, beside the others:

```fortran
use mod_test_tox_foo, only: get_all_tests_tox_foo
...
call add_suite("tox_foo", get_all_tests_tox_foo)
```

To add a case, write its subroutine, raise the `allocate` count and add its `test_case` line.

### A suite with one module per procedure

A module with several published procedures gets a directory, `test/mod_test_<module>/`, holding:
- the suite module `mod_test_<module>.F90`, which only concatenates its children's arrays (see
  `mod_test_tox_normalization.F90`);
- one child per procedure, `mod_test_<procedure>.F90`, each with its own
  `get_all_tests_<procedure>()`. The procedure name alone is enough: published names are unique
  across the library, because Python's namespace is flat, and the directory already names the
  module. With the module in front, `tox_relative_axis_plane_tools`' children would exceed
  Fortran's 63-character limit for names;
- data the children share in a module of their own (e.g. `mod_test_tox_normalization_fixtures.F90`).

`run_tests` registers only the suite module; it never sees the children. Test-case names must be
unique within the suite.

### What makes a good case

- **Hand-derived expected values.** An expectation computed with the code's own formula cannot
  fail. Prefer inputs whose results are exact, and say in the comment how you got them.
- **Every early exit, bound and dimension on its own.** Setting all dimensions to zero at once
  exercises only the first check. Test both sides of each bound: the last valid value and the
  first invalid one.
- **The documented errors.** Compare `get_err_code(ierr)` with an `ERR_*` name, never a literal.
  `assert_err(ierr, ERR_X, msg, arg_pos)` also checks which argument was blamed.
- **NaN and Inf are rejected.** TOX takes no NaN or Inf and writes none, except in preprocessing
  such as data integration. A procedure's suite checks that it rejects both.
- **No array temporaries.** Assign array constructors to variables rather than passing
  expressions; `--diagnostics` warns about temporaries.
- **Assertions that can fail.** A size fixed at compile time, or `no NaN` after an exact
  comparison, checks nothing.

### Assertions

All in `test_framework/src/asserts.F90`; each takes the failure message as its last argument (before
optional ones).

| Assertion | Checks |
|---|---|
| `assert_true(cond, msg)`, `assert_false(cond, msg)` | a condition, default or `c_bool` logical |
| `assert_equal_int(a, b, msg)`, `assert_not_equal_int` | integers |
| `assert_equal_real(actual, expected, tol, msg)`, `assert_not_equal_real` | reals, within an absolute tolerance |
| `assert_equal_array_real(actual, expected, n, tol, msg)` | the first `n` elements; any rank, by sequence association |
| `assert_allclose_array_real` | arrays within a tolerance |
| `assert_equal_array_int`, `assert_equal_array_char`, `assert_equal_array_logical` | arrays of other types |
| `assert_equal_complex`, `assert_not_equal_complex`, `assert_equal_array_complex` | complex values |
| `assert_in_range_real(value, min, max, msg)`, `assert_in_range_int` | a value within bounds |
| `assert_no_nan_real(array, n, msg)`, `assert_no_inf_real` | no NaN / no Inf |
| `assert_sorted_real`, `assert_sorted_int`, `assert_unique_int`, `assert_permutation` | ordering properties |
| `assert_contains_int`, `assert_array_int_contains`, `assert_sum_equal` | membership and sums |
| `assert_string_equal`, `assert_string_contains` | strings |
| `assert_err(ierr, expected_code, msg, arg_pos)` | an error code, and optionally the argument it blames |

## Writing Python and R tests

A binding suite is `python/test/mod_test_<module>.py` or `r/test/mod_test_<module>.R`. Every
function named `test_*` in it is a case:

```python
def test_calc_fchange():
    result = calc_fchange(np.array([1, 1], dtype=np.int32), np.array([2, 3], dtype=np.int32), averages)
    _assert_matrix(result, (2, N_GENES), "calc_fchange")

def test_calc_fchange_rejects_a_tissue_past_the_last():
    assert_error(lambda: calc_fchange(...), "there are only three tissues", ERR_INVALID_INPUT)

if __name__ == '__main__':
    run_all_tests(globals().values())
```

In R the file ends with `run_all_tests()`, and the error check is
`assert_error(expr, msg, ERR_INVALID_INPUT)`. Omitting the code demands an error that is not a
TensorOmics error, such as the binding's own type or shape check. Run a suite from the repository
root, `python3 python/test/mod_test_tox_normalization.py` or
`Rscript r/test/mod_test_tox_normalization.R`, after a build; `run_all_tests.sh` runs them all.
