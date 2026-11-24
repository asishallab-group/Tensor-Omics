from pathlib import Path

import pytest

from codegen.config import CONVENTIONS, Conventions, Paths


@pytest.mark.parametrize(
    "name, expected",
    [
        ("mode", "mode"),
        ("method", "method"),
        ("MODE", "mode"),
        ("link_method", "method"),
        ("clustering_mode", "mode"),
        # a name merely containing an alias is not a mode argument
        ("modern", None),
        ("methodology", None),
        ("node", None),
        # the alias must be a whole trailing word
        ("linkmethod", None),
        ("n_modes", None),
        ("vector", None),
        # a bare suffix has no owning name in front of it
        ("_mode", None),
    ],
)
def test_mode_alias_of(name, expected):
    assert CONVENTIONS.mode_alias_of(name) == expected


@pytest.mark.parametrize(
    "name, expected",
    [
        ("genes_mask", "genes"),
        ("genes_selection_mask", "genes"),
        ("GENES_MASK", "genes"),
        # the longer suffix wins, so the owner is 'genes' and not 'genes_selection'
        ("genes_selection_mask", "genes"),
        # a bare suffix has no owner to belong to
        ("_mask", None),
        ("mask", None),
        ("_selection_mask", None),
        ("masked_genes", None),
        ("vector", None),
    ],
)
def test_mask_arg_name_of(name, expected):
    assert CONVENTIONS.mask_arg_name_of(name) == expected


def test_mask_suffixes_are_ordered_longest_first():
    # otherwise '_mask' would match 'genes_selection_mask' and claim the wrong owner
    suffixes = CONVENTIONS.mask_suffixes
    assert list(suffixes) == sorted(suffixes, key=len, reverse=True)


def test_conventions_are_immutable():
    with pytest.raises(Exception):
        CONVENTIONS.temporary_prefix = "temp_"


def test_conventions_can_be_varied_for_tests():
    variant = Conventions(temporary_prefix="work_")

    assert variant.temporary_prefix == "work_"
    assert CONVENTIONS.temporary_prefix == "tmp_"


def test_paths_resolve_relative_against_root():
    paths = Paths(root=Path("/repo"))

    assert paths.resolve(paths.src_dir) == Path("/repo/src")


def test_paths_leave_absolute_paths_alone():
    paths = Paths(root=Path("/repo"))

    assert paths.resolve(Path("/elsewhere/src")) == Path("/elsewhere/src")


def test_src_dir_is_injectable_so_tests_can_parse_fixtures():
    paths = Paths(root=Path("/repo"), src_dir=Path("fixtures"))

    assert paths.resolve(paths.src_dir) == Path("/repo/fixtures")


@pytest.mark.parametrize(
    "name, expected",
    [
        ("bitmask", True),
        ("bitmasks", True),
        ("family_bitmasks", True),
        ("filtered_paralogs_bitmask", True),
        ("Source_BitMask", True),
        # the word must stand alone or follow an underscore
        ("tmp_bitmasker", False),
        ("nobitmask", False),
        # the selection-mask suffixes are another role
        ("genes_mask", False),
        ("bit_mask", False),
    ],
)
def test_is_bitmask_name(name, expected):
    assert CONVENTIONS.is_bitmask_name(name) is expected


@pytest.mark.parametrize("name", ["bitmask", "family_bitmasks", "filtered_paralogs_bitmask"])
def test_a_bitmask_is_never_a_selection_mask(name):
    # `_bitmask` does not end in `_mask`, so `n_selected_<x>` can never count one
    assert CONVENTIONS.mask_arg_name_of(name) is None
