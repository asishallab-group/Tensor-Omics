#!/usr/bin/env python3
"""The package's top level exports the error module whole (#221 WP7): every ERR_* and STAT_*
code, ARG_POS_FACTOR, every exception class and check_err_code, so a caller compares a code or
catches by class without knowing the submodule."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import tensor_omics
from tensor_omics import (
    ARG_POS_FACTOR,
    ERR_DIVISION_BY_ZERO,
    ERR_OK,
    ToxArgumentError,
    normalize_all_trajectories,
)
from test_helpers import run_all_tests


def test_the_status_decodes_with_the_names_the_package_exports():
    """A status array holds tox codes; the package itself names them"""
    trajectories = np.ones((3, 4, 5), dtype=np.float64) * 7.0
    trajectories[1, :, :] = np.arange(5).reshape(1, 1, 5)  # only factor 1 varies

    status = normalize_all_trajectories(trajectories)["status"]

    codes = status % ARG_POS_FACTOR
    assert np.all(codes[[0, 2], :] == ERR_DIVISION_BY_ZERO), "a constant factor cannot be scaled"
    assert np.all(codes[1, :] == ERR_OK), "a varying factor normalizes"


def test_the_binding_s_own_refusal_is_an_exported_class():
    try:
        normalize_all_trajectories("not an array")
    except ToxArgumentError:
        pass
    else:
        raise AssertionError("a string for an array should be refused by the binding")


def test_every_exported_name_is_the_error_module_s_own():
    from tensor_omics import error_handling
    for name in ("ERR_DIVISION_BY_ZERO", "ARG_POS_FACTOR", "ToxArgumentError", "check_err_code"):
        assert name in tensor_omics.__all__, f"{name} is not in __all__"
        assert getattr(tensor_omics, name) is getattr(error_handling, name), name


if __name__ == "__main__":
    run_all_tests(globals().values())
