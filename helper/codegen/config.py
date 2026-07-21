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
    c_interface_category: str = "C-interface"

    #: Suffix of the allocating variant of a procedure pair
    alloc_suffix: str = "_alloc"
    #: Infix for the non-allocating variant, used only when an `_alloc` sibling exists
    expert_infix: str = "_expert"
    #: Suffix of every generated C symbol and of the modules holding them
    c_suffix: str = "_c"

    #: Prefix marking a work array: allocated silently, never returned
    temporary_prefix: str = "tmp_"
    #: Suffix marking the argument that carries another argument's shape
    shape_suffix: str = "_shape"
    #: Prefix of the count belonging to a `<arg>_mask` / `<arg>_selection_mask`
    mask_count_prefix: str = "n_selected_"
    #: Suffixes recognised on a mask argument, longest first so matching is unambiguous
    mask_suffixes: tuple[str, ...] = ("_selection_mask", "_mask")

    #: Accepted spellings of a mode argument. A dummy is a mode argument if its name is
    #: an alias verbatim (`mode`) or carries one as a suffix (`link_method`).
    mode_aliases: tuple[str, ...] = ("mode", "method")

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
    c_interface_dir: Path = Path("src/c_interface")
    python_out_dir: Path = Path("python/tensor_omics")
    rcpp_out_dir: Path = Path("rcpp/tensor_omics")

    def resolve(self, path: Path) -> Path:
        return path if path.is_absolute() else self.root / path


CONVENTIONS = Conventions()
