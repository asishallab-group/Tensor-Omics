# Code generator footgun: `n_selected_<X>` silently loses its mask if the mask isn't named `<X>_mask`/`<X>_selection_mask`

## Where this lives in the generator

`helper/codegen/ir/roles.py`, `_link_mask_counts`/`_find_mask_for`, together with
`helper/codegen/config.py`'s `Conventions.mask_arg_name_of`:

```python
# roles.py
def _link_mask_counts(self) -> None:
    prefix = self.conventions.mask_count_prefix          # "n_selected_"
    for argument in self.arguments:
        name = argument.name.lower()
        if not name.startswith(prefix):
            continue
        owner_name = name[len(prefix):]                  # "n_selected_seed" -> "seed"
        if not owner_name:
            continue
        mask = self._find_mask_for(owner_name)
        if mask is None:
            continue                                      # <-- silently gives up here
        argument.roles.mask_count_of = mask
        mask.roles.count_arg = argument

def _find_mask_for(self, owner_name: str) -> Argument | None:
    for candidate in self.arguments:
        if self.conventions.mask_arg_name_of(candidate.name) == owner_name:
            return candidate
    return None

# config.py
mask_suffixes: tuple[str, ...] = ("_selection_mask", "_mask")

def mask_arg_name_of(self, name: str) -> str | None:
    lowered = name.lower()
    for suffix in self.mask_suffixes:
        if lowered.endswith(suffix):
            owner = lowered[: -len(suffix)]
            if _is_identifier(owner):
                return owner
    return None
```

The link is found by **exact string equality** between the count's derived owner
(`n_selected_<owner>` -> `<owner>`) and the mask's derived owner (`<owner><suffix>` ->
`<owner>`, suffix-stripped). If the two don't produce the identical `<owner>` string, no link
is made -- and nothing in `_link_mask_counts` reports that.

## The problem

Section 5.3 of `codegen_guide.md` documents the convention as `<arg>_mask` /
`n_selected_<arg>`, with a worked example (`axes_selection_mask` / `n_selected_axes`). It
reads as a naming *pattern* to follow loosely. In practice the match is a **literal suffix
strip and string-equality check**, not a semantic association -- so any mask name that isn't
*exactly* `<owner>_mask` or `<owner>_selection_mask` fails to link, even when a human reading
the two argument names side by side would consider the pairing obvious.

## Minimal example

```fortran
pure subroutine gene_stats_kernel(gene_values, n_genes, is_gene_mask, n_selected_gene, result)
    real(real64),   intent(in)  :: gene_values(n_genes)
    integer(int32), intent(in)  :: n_genes
    logical,        intent(in)  :: is_gene_mask(n_genes)
        !! .TRUE. for genes to include
    integer(int32), intent(in)  :: n_selected_gene
        !! Number of selected genes (count of .TRUE. in is_gene_mask)
    real(real64),   intent(out) :: result(n_selected_gene)
    ...
end subroutine gene_stats_kernel
```

`n_selected_gene` strips to owner `"gene"`. `is_gene_mask` strips `_mask` to owner
`"is_gene"`, not `"gene"` -- the leading `is_` is included in what gets compared, since
`mask_arg_name_of` only strips the trailing `_mask`/`_selection_mask` suffix, nothing at the
front. `"is_gene" != "gene"`, so `_find_mask_for("gene")` returns `None`, and
`_link_mask_counts` just `continue`s.

Renaming the mask to `gene_mask` (or `gene_selection_mask`) -- nothing else -- restores the
link.

## Effect of the silent miss

Two things regress at once, neither reported anywhere:

1. **The cross-check the whole convention exists for is skipped.** With the link made, the
   generated wrapper emits `if (count(gene_mask, kind=int32) /= n_selected_gene)
   call set_err_once(...)`. Without it, that line is never generated -- a caller can pass a
   `n_selected_gene` that doesn't match `count(is_gene_mask)` at all, and the wrapper will not
   notice.
2. **`n_selected_gene` falls back to generic extent validation** (`validate_dimension_size`,
   which rejects 0 as `ERR_EMPTY_INPUT` and otherwise just checks non-negativity) instead of
   whatever `DM_MIN`/`DM_MAX` the argument itself documents, unless those are given
   explicitly. In our case (a kernel that legitimately wants `n_selected_seed = 0` to be a
   valid, empty-result input) this cost us a second silent surprise.

For a hand-written kernel that then does its own unchecked bookkeeping keyed on
`n_selected_<X>` -- e.g. sizing a local array to `n_selected_gene` and filling it by scanning
the mask -- a caller who *does* manage to pass a mismatched count (trivial for a Fortran
caller using the expert/kernel tier directly, which skips the validating wrapper entirely by
design) gets undefined behavior: an out-of-bounds write if the true count is larger than
`n_selected_gene`, or uninitialized trailing entries if it's smaller. The generator gave no
signal at generation time (`--check` reports 0 warnings either way) that the safety net meant
to prevent exactly this was never wired up.

## Suggested fixes (either would have caught this)

1. **Warn when a count has no mask.** In `_link_mask_counts`, when `mask is None`, emit a
   diagnostic (warning, or error under `--check`'s "0 warnings" bar) naming the expected mask
   spelling(s): `"'n_selected_gene' looks like a mask count, but no 'gene_mask' or
   'gene_selection_mask' argument was found"`. This is the cheap fix -- it doesn't change
   matching behavior at all, it just stops the miss from being silent. `_link_result_sizes`,
   a few lines down in the same file, already does exactly this pattern (`diagnostics.error`
   when a `DM_RESULT_SIZE_IS` target can't be resolved) -- `_link_mask_counts` is the odd one
   out for staying silent.
2. **Make the pairing explicit and generator-checked**, the same way `DM_OUTPUT_FROM` makes a
   workspace-size dependency explicit rather than inferred: something like
   `DM_MASK_FOR(is_gene_mask)` on the count argument, or `DM_COUNT_OF(n_selected_gene)` on the
   mask argument, checked against the procedure's actual argument list at generation time (an
   unresolvable reference is already an error elsewhere in the generator, e.g.
   `DM_RESULT_SIZE_IS`'s own diagnostic above). This removes the naming coincidence
   requirement entirely, at the cost of one more annotation to write.

Either is compatible with keeping the existing bare `<arg>_mask` / `n_selected_<arg>`
convention as the zero-annotation default for the common case -- the fix is only about making
a miss visible, not about changing what matches today.

---

# Second footgun: a literal `/` in an array-bound expression breaks the Python binding (integer division vs. true division)

## Where this lives in the generator

`helper/codegen/emit/python_ctypes.py`, `_python_extent`:

```python
_FORTRAN_KIND_SUFFIX = re.compile(r"\b(\d+(?:\.\d+)?)_[A-Za-z]\w*")
_FORTRAN_SIZE = re.compile(r"\bsize\s*\(\s*([A-Za-z_]\w*)\s*(?:,\s*(\d+)\s*)?\)")

def _python_extent(extent: str) -> str:
    """A Fortran extent expression as Python.

    `max`/`min` are built-ins already. What needs translating is the numeric kind suffix,
    which Python has no notion of, and `size(...)`, which numpy spells as an attribute.
    """
    return _FORTRAN_SIZE.sub(_size_as_numpy, _FORTRAN_KIND_SUFFIX.sub(r"\1", extent))
```

used in `_new_value` to build the shape tuple for `np.empty`:

```python
shape = ", ".join(_python_extent(e) for e in argument.dimension.extents)
return f"np.empty(({shape},), dtype={dtype_of(argument)}, order={order})"
```

## The problem

`_python_extent` only rewrites two things: the `_int32`-style kind suffix, and `size(...)`
calls. Every operator in the expression -- `+`, `-`, `*`, `/`, `min`, `max` -- passes through
completely unchanged, on the assumption that Fortran and Python spell arithmetic the same way.
They don't, for one operator: Fortran's `/` between two `integer` operands is truncating
integer division; Python's `/` is always true (float) division, unconditionally, regardless of
operand type. `//` is Python's integer-division spelling, and nothing produces it here.

## Minimal example

```fortran
integer(int32), intent(in) :: n
    !! DM_MIN(2_int32)
integer(int32), intent(out) :: pairs(n*(n - 1)/2)
    !! One row per unordered pair of the n items
```

In Fortran, `n*(n - 1)/2` is exact integer arithmetic (the product of two consecutive
integers is always even, so the truncating division loses nothing) -- this is the ordinary way
to spell a triangular number as an array bound. The generated Python binding computes the
same shape expression almost verbatim:

```python
pairs = np.empty((n*(n-1)/2,), dtype=np.int32, order='F')
```

`n*(n-1)/2` is a `float` in Python (e.g. `6*5/2` is `15.0`, not `15`). `np.empty` requires an
integer shape, so this raises `TypeError: 'float' object cannot be interpreted as an integer`
-- for every call, unconditionally, regardless of what `n` actually is. `--check` reports 0
warnings; the break only surfaces the first time the Python binding is actually exercised.

## Where we hit this

`ensemble_reconciliation_kernel` in
`src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_reconciliation_kernel.F90`
originally sized its output at `n_ensembles*(n_ensembles - 1)/2`, the exact worst-case count
of unordered ensemble pairs. The Fortran kernel, the generated Fortran wrapper, and the R
binding (which never recomputes the shape itself -- it reads the buffer `.Call` already
allocated on the Fortran/C side) all worked fine. Only the Python binding's own
pre-allocation broke. Worked around by doubling the bound to the division-free
`n_ensembles*(n_ensembles - 1)` (a safe, if unnecessarily loose, upper bound) instead --
correct in all three languages, at the cost of an avoidable 2x over-allocation.

## Suggested fixes

1. **Have `_python_extent` (and its R counterpart, if the R backend ever starts
   recomputing shapes instead of trusting the Fortran side) rewrite integer division.**
   Fortran's own type-checking already knows whether both operands of a `/` are integer typed
   at the point the extent expression is parsed -- if that type information is available where
   `_python_extent` runs, rewrite `/` to `//` (and the equivalent `%/%` in R) exactly when both
   sides are integer; leave it alone otherwise (a `real` bound is legitimate and division
   there does mean true division in both languages).
2. **Failing that, warn.** If plumbing per-operand type information into `_python_extent` is
   more invasive than it's worth, at minimum scan a kernel's array-bound expressions for a
   bare `/` between what look like integer-typed operands and flag it under `--check` --
   `"'pairs' is sized 'n*(n - 1)/2': Fortran integer division here becomes Python float
   division in the generated binding, which np.empty will reject at every call. Use `//`
   -- oh wait, this is Fortran, not Python -- rephrase without a literal division, e.g. by
   doubling the bound"`. Even a warning that just says "a `/` here will not round-trip to
   Python the way you expect" saves the same trial-and-error this took to track down.
3. **Document the gap explicitly** in `codegen_guide.md` section 5 wherever array-bound
   specification expressions are discussed, so a kernel author knows to avoid a bare integer
   `/` in a bound expression (or knows why doubling it, as done here, is a legitimate
   workaround) without having to hit the `TypeError` first.

---

# Third footgun: an omitted optional argument skips validation entirely, even when its `DM_DEFAULT` can violate a runtime-dependent `DM_MAX`

## Where this lives in the generator

The generated wrapper's own validation call, e.g.
`src/generated/tox/tox_shape_truthful_clustering_seeding.F90`:

```fortran
call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
```

and `validate_in_range_int` itself, `src/f42/tox_errors.F90`:

```fortran
pure subroutine validate_in_range_int(val, ierr, arg_pos, min, max, sentinel)
    integer(int32), intent(in), optional :: val
    ...
    if (present(val)) then
        ...
        if ((val < actual_min) .or. (val > actual_max)) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
        end if
    end if
end subroutine validate_in_range_int
```

## The problem

The wrapper passes its own (still-`optional`) `k_density` dummy straight through to
`validate_in_range_int`'s own `val`, which is *also* `optional` and gates the entire check
behind `if (present(val))`. If the original caller omitted `k_density`, it is never present
at that call site either, so the bound check is skipped -- not "skipped because it's fine,"
skipped as in *never evaluated at all*. The unresolved-absent optional is then passed on to
the kernel unchanged, where `M_DEFAULT_VAL` resolves it to the literal `DM_DEFAULT` constant
(e.g. 30) with no bound check ever having run.

This is fine whenever the default constant is a legitimate value for every possible input --
but `DM_MAX(n_vectors - 1_int32)` is not a constant, it depends on the *other* arguments of
this specific call. A `DM_DEFAULT` of 30 is not a valid value once `n_vectors <= 30`, and
nothing catches that.

## Minimal example

```fortran
pure subroutine density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                      k_density, tmp_neighbors, tmp_distances, tmp_range_stack, ..., labels)
    ...
    integer(int32), intent(in), optional :: k_density
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors - 1_int32)
        !! DM_DEFAULT(30_int32)
    integer(int32), intent(out) :: tmp_neighbors(n_vectors)
    ...
    integer(int32) :: actual_k_density, k_query
    M_DEFAULT_VAL(k_density, actual_k_density, 30_int32)
    k_query = actual_k_density + 1
    call kd_knn_query_helper(..., tmp_neighbors(1:k_query), ...)   ! n_vectors=5, k_query=31
```

Calling the *generated, validating* wrapper with `k_density` omitted, on a 5-point dataset,
never raises `ERR_INVALID_INPUT` -- it resolves `k_density` to 30 inside the kernel and asks
`kd_knn_query_helper` to fill a 31-element slice of a 5-element array. That is an
out-of-bounds write. We hit this directly: it corrupted memory silently, then crashed later
in unrelated code (a `SIGSEGV` inside `libgfortran`'s `_dl_fini`, during process teardown,
nowhere near the actual bug) -- the kind of failure that is very expensive to trace back to
its real cause.

## Why this is easy to miss

Every *explicit* value is validated correctly -- `test_growth_radius_k_min_too_large` and
`test_seeds_invalid_k_density` (both passing an explicit, too-large value) both correctly get
`ERR_INVALID_INPUT`. The gap is specifically the combination of (a) omitting the optional and
(b) a dataset smaller than the default. Python and R callers are not exposed to this in
practice: both bindings always resolve and pass an explicit value at the C-ABI boundary
(confirmed by reproducing the equivalent call from Python, which raised a clean
`ToxInputError` rather than crashing) -- a Python/R caller never actually presents Fortran
with a truly-absent optional. The risk is specifically a *Fortran* caller (another kernel, a
hand-written test, a future orchestrator) that omits the argument on data smaller than the
default -- exactly what one of our own first-draft tests did by accident.

## Suggested fixes

1. **Resolve the default before validating, not after.** If the wrapper computed
   `actual_k_density` via the same `M_DEFAULT_VAL` logic the kernel uses, then validated
   *that* (always present) against `DM_MIN`/`DM_MAX`, a default that violates a
   runtime-dependent bound would be caught for every call, not just the ones that happen to
   pass an explicit value. This also matches what a reader would reasonably assume
   `DM_DEFAULT` already means: "a value that passes validation," not "a value substituted
   after validation had its chance."
2. **Failing that, warn at generation time when this is even possible.** If a `DM_MAX` (or
   `DM_MIN`) expression references another argument (not a compile-time constant) *and* the
   argument carries a `DM_DEFAULT`, that combination can never be verified as safe at
   generation time -- flag it under `--check`, so a kernel author at least knows the default
   is unchecked at the point they write it, rather than discovering it via a crash in
   unrelated code.

Neither fix has been made in the generator itself yet -- both are still open. What *has* been
fixed, at the two call sites this actually bit us
(`calc_ensemble_growth_radius_kernel`'s `k_min`, `density_labels_kernel`'s `k_density`), is a
kernel-level stopgap: `actual_k = min(actual_k, n_vectors - 1)` right after `M_DEFAULT_VAL`,
which makes the *resolved* value safe regardless of whether it came from an explicit
(already wrapper-validated) argument or the unchecked default. This is a per-kernel patch,
not a fix to the generator -- any *other* kernel with the same `DM_DEFAULT` +
runtime-dependent-`DM_MAX` combination remains exposed until either it gets the same local
clamp or the generator itself is fixed per the two suggestions above. Regression tests for
the two fixed call sites:
`test/mod_test_shape_truthful_clustering_seeding.F90::test_density_labels_omitted_k_density_is_clamped`/
`test_seeds_omitted_k_density_is_clamped`, and
`test/mod_test_shape_truthful_clustering_ensemble_growing.F90::test_growth_radius_omitted_k_min_is_clamped`.
