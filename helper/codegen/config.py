"""Central configuration for the code generator.

Everything the generator recognises by convention lives here, so a convention is
changed in one place instead of being grepped for across the emitters.

Two kinds of settings are distinguished:

- `Conventions`: the source-language contract (prefixes, suffixes, meta tags). These
  are what a Fortran author writes and what `README.md` documents. Changing them is an
  API change for the whole codebase, so they are not meant to be overridden per run --
  they are grouped as a class purely so tests can construct a variant.
- `Paths`: where the generator reads and writes. Injectable, because the test suite
  points the frontend at fixture sources instead of `src/`.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

#: A Fortran identifier: a letter followed by letters, digits or underscores.
_IDENTIFIER_RE = re.compile(r"[a-z][a-z0-9_]*\Z", re.IGNORECASE)


def _is_identifier(name: str) -> bool:
    return _IDENTIFIER_RE.match(name) is not None


@dataclass(frozen=True)
class Conventions:
    """The naming/annotation contract that Fortran sources must follow."""

    #: The Ford `category` value marking a procedure for export. This is a fallback for a
    #: hand-built Conventions; a real run reads it from the `M_EXPORT_C` macro (see
    #: `frontend.export_category`), so the marker and the generator cannot drift.
    c_binding_category: str = "C-binding"

    #: Suffix of the allocating variant of a procedure pair
    alloc_suffix: str = "_alloc"
    #: Infix for the non-allocating variant, used only when an `_alloc` sibling exists
    expert_infix: str = "_expert"
    #: Suffix of every generated C symbol and of the modules holding them
    c_suffix: str = "_c"

    #: Suffix marking a hand-written kernel procedure (and its module). The generator
    #: turns `loess_fit_kernel` in module `tox_loess_kernel` into the validating wrapper
    #: `loess_fit` (and, when it needs work arrays, `loess_fit_alloc`) in module `tox_loess`.
    kernel_suffix: str = "_kernel"
    #: Suffix of the third variant, the plain validating wrapper. Empty -- it takes the
    #: bare name -- but named so the three variants (kernel, validating, alloc) are
    #: symmetric and the naming lives in one place.
    validating_suffix: str = ""
    #: Suffix marking a permutation vector. In the allocating wrapper a `<base>_perm`
    #: argument is seeded with `init_perm` and heapsorted against `<base>`.
    perm_suffix: str = "_perm"

    #: Prefix marking a work array: allocated silently, never returned
    temporary_prefix: str = "tmp_"
    #: Suffix marking the argument that carries another argument's shape
    shape_suffix: str = "_shape"
    #: Prefix of the count belonging to a `<arg>_mask` / `<arg>_selection_mask`
    mask_count_prefix: str = "n_selected_"
    #: Suffixes recognised on a mask argument, longest first so matching is unambiguous
    mask_suffixes: tuple[str, ...] = ("_selection_mask", "_mask")
    #: A square real matrix whose name is one of these, or ends in `_<one>`, is validated for
    #: distance-matrix structure (symmetry, non-negativity, zero diagonal) rather than finiteness
    distance_matrix_suffixes: tuple[str, ...] = ("distances", "distance_matrix")

    #: Accepted spellings of a mode argument. A dummy is a mode argument if its name is
    #: an alias verbatim (`mode`) or carries one as a suffix (`link_method`).
    mode_aliases: tuple[str, ...] = ("mode", "method")

    #: The flag a prologue sets to say it handled the call, so the kernel is skipped
    prologue_handled_arg: str = "handled"

    #: Header of the optional third mode-table column. Its presence opts a mode argument
    #: into per-mode wrapper generation: one procedure per mode, named in that column
    #: (e.g. `detect_dosage_effect`), instead of a single procedure taking the mode.
    mode_procedure_header: str = "procedure"

    #: Header of the first column of the table that maps a DM_OUTPUT_FROM producer's
    #: input names onto the consumer's arguments, where the two differ
    producer_input_header: str = "producer input"
    #: Header of its second column
    producer_supplied_by_header: str = "supplied by"

    #: Prefix of an error code parameter in `tox_errors`; these raise in Python/R
    error_code_prefix: str = "ERR_"
    #: Prefix of a status code parameter; these are outcomes, not failures, and never raise
    status_code_prefix: str = "STAT_"
    #: The error code parameter meaning "no error"
    ok_code: str = "ERR_OK"

    #: Name of the error argument every C wrapper carries, synthesised when absent
    error_arg: str = "ierr"

    def mode_alias_of(self, name: str) -> str | None:
        """Return the mode alias `name` uses, or None if it is not a mode argument.

        Matches `mode`, `method`, and any `*_mode` / `*_method` spelling.
        """
        lowered = name.lower()
        for alias in self.mode_aliases:
            if lowered == alias:
                return alias
            if lowered.endswith(f"_{alias}") and _is_identifier(lowered[: -len(alias) - 1]):
                return alias
        return None

    def mask_arg_name_of(self, name: str) -> str | None:
        """Return the argument a mask belongs to, e.g. `foo_selection_mask` -> `foo`.

        The remainder has to be a Fortran identifier in its own right, otherwise
        `_selection_mask` would match the shorter `_mask` suffix and name a
        non-existent `_selection` argument as its owner.
        """
        lowered = name.lower()
        for suffix in self.mask_suffixes:
            if lowered.endswith(suffix):
                owner = lowered[: -len(suffix)]
                if _is_identifier(owner):
                    return owner
        return None


@dataclass(frozen=True)
class Paths:
    """Filesystem layout. Relative paths are resolved against `root`."""

    root: Path = field(default_factory=Path)
    src_dir: Path = Path("src")
    macros_header: Path = Path("src/macros.h")
    c_binding_dir: Path = Path("src/bindings/c")
    r_binding_dir: Path = Path("src/bindings/r")
    #: Hand-written kernels; the generated wrappers derived from them land in `tox_out_dir`
    kernel_src_dir: Path = Path("src/kernel")
    #: Where the generated `foo` / `foo_alloc` wrapper modules are written (generated-only)
    tox_out_dir: Path = Path("src/tox")
    python_out_dir: Path = Path("python/tensor_omics")
    r_out_dir: Path = Path("r/tensor_omics")
    snippets_dir: Path = Path("snippets")

    def resolve(self, path: Path) -> Path:
        return path if path.is_absolute() else self.root / path


CONVENTIONS = Conventions()
