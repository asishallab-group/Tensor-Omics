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
    #: The doc macro an author writes to mark a procedure for export. The category above is
    #: read from it, so this is the name the diagnostics quote back.
    export_marker: str = "M_EXPORT_C"

    #: Suffix of the allocating half of a *hand-written* procedure pair. The generator never
    #: emits it: a generated allocating wrapper takes the plain name. It is still recognised
    #: because f42 writes its pairs by hand (`compute_edf` / `compute_edf_alloc`), and that is
    #: the shape `stripped_name` translates into the published `foo` / `foo_expert`.
    alloc_suffix: str = "_alloc"
    #: Suffix of the non-allocating entry point, in Fortran and as the binding languages
    #: publish it. A generated pair is named `foo` / `foo_expert` outright; a hand-written
    #: `foo` / `foo_alloc` pair is published under those names by `abi.c_abi.stripped_name`.
    expert_suffix: str = "_expert"
    #: Suffix of every generated C symbol and of the modules holding them
    c_suffix: str = "_c"

    #: Suffix marking a hand-written implementation procedure (and its module). The generator
    #: turns `loess_fit_impl` in module `tox_loess_impl` into the entry point `loess_fit` in
    #: module `tox_loess` -- and, where it has work arrays to take over, into the pair
    #: `loess_fit` (which allocates) and `loess_fit_expert` (which does not).
    impl_suffix: str = "_impl"
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

    #: The flag a prologue sets to say it handled the call, so the implementation is skipped
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
    #: Everything this generator writes lives under `src/generated`, and nothing hand-written
    #: does. That is the whole rule a reader (or a review, or an editor's ignore list) needs.
    #: The wrapper modules mirror the tree they came from into it: an implementation at
    #: `src/<rest>/<module>_impl.F90` generates `src/generated/<rest>/<module>.F90`, so the
    #: layer a procedure belongs to survives generation and no directory is special-cased.
    generated_dir: Path = Path("src/generated")
    c_binding_dir: Path = Path("src/generated/bindings/c")
    r_binding_dir: Path = Path("src/generated/bindings/r")
    python_out_dir: Path = Path("python/tensor_omics")
    r_out_dir: Path = Path("r/tensor_omics")
    snippets_dir: Path = Path("snippets")

    def resolve(self, path: Path) -> Path:
        return path if path.is_absolute() else self.root / path


CONVENTIONS = Conventions()
