r"""tox_data_integration_js_comp_test_permutation_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test

Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
(without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
each study's observed global JSD with the values the shuffled data sets give.

**Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
numbered in slot order, replicate fastest, then study.

**Entry points.** One assignment of samples to slots, one JS-Comp-Test:
:func:`tensor_omics.run_js_comp_test_on_sample_assignment`
(k-NN neighborhoods, the fixed `(n_points, n_neighbors)` of `run_js_comp_test`) and
:func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
(adaptive neighborhoods, `run_js_comp_test_adaptive`). The whole test, `n_permutations` random
assignments plus the observed one:
:func:`tensor_omics.run_js_comp_test_label_permutation_test`
and
:func:`tensor_omics.run_js_comp_test_adaptive_label_permutation_test`.
`calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds` size
the work arrays for a caller that allocates its own.

**What a p-value means.** Study `i`'s p-value is `(count + 1)/(n_valid + 1)`, where `count` is
how many of the `n_valid` valid shuffled data sets give study `i` a global JSD at least as large
as the observed one (the +1 counts the observed data itself, so it is never 0). A small p-value
(e.g. 0.01) says: hardly any random relabelling makes study `i` look as different from the
others as the real labels do, so its samples differ systematically (a batch effect). A large
one says the real labels are nothing special.

This test is not the residual-level permutation test of `run_js_comp_test`
(:func:`tensor_omics.gjct_permutation_test`), which keeps
the data and the histograms fixed and only reshuffles the pooled histogram counts.

Python binding, generated from tox_data_integration_js_comp_test_permutation_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_js_comp_test_permutation_bounds_c.restype = None
_lib.calc_js_comp_test_permutation_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_points", "n_neighbors", "n_bootstraps", "two_sided_bootstrapping_significance_level", "max_n_points_candidate", "max_n_neighbors_candidate", "n_bootstrapping_top_k_jsds", "ierr",)

_lib.calc_adaptive_js_comp_test_permutation_bounds_c.restype = None
_lib.calc_adaptive_js_comp_test_permutation_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ADAPTIVE_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "k_start", "n_bootstraps", "two_sided_bootstrapping_significance_level", "max_n_points_candidate", "n_bootstrapping_top_k_jsds", "ierr",)

def calc_js_comp_test_permutation_bounds(
        max_n_genes_all_studies,
        n_points,
        n_neighbors,
        n_bootstraps=1000,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the work-array capacities of the k-NN study-label Permutation-Test

    Sizes the per-point and per-neighbor work arrays of
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`
    and
    :func:`tensor_omics.run_js_comp_test_label_permutation_test`
    so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
    `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
    K)`, with `P` and `K` the bounds of
    :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
    `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    :func:`tensor_omics.calc_js_comp_test_n_top_k_jsds` computes it.

    Rejects an out-of-range argument with invalid input (a non-finite significance level with
    `ERR_NAN_INF`); the outputs are then left undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_points : int
        Number of reference points of the fixed-parameter run
        The minimum valid value is `1`.
    n_neighbors : int
        Number of neighbors per neighborhood of the fixed-parameter run
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps of the parameter search, when it is rerun
        The minimum valid value is `1`.
        The default value is `1000`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Capacity for reference points
        max_n_neighbors_candidate : int
            Capacity for neighbors per neighborhood
        n_bootstrapping_top_k_jsds : int
            Bootstrap heap size (values kept at each end of the bootstrap distribution)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_impl::calc_js_comp_test_permutation_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    max_n_neighbors_candidate = ctypes.c_int(0)
    n_bootstrapping_top_k_jsds = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_permutation_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(max_n_neighbors_candidate),
        ctypes.byref(n_bootstrapping_top_k_jsds),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
        "n_bootstrapping_top_k_jsds": n_bootstrapping_top_k_jsds.value,
    }

def calc_adaptive_js_comp_test_permutation_bounds(
        max_n_genes_all_studies,
        n_studies,
        k_start,
        n_bootstraps=1000,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the work-array capacity of the adaptive study-label Permutation-Test

    Sizes the per-point work arrays of
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    and
    :func:`tensor_omics.run_js_comp_test_adaptive_label_permutation_test`:
    `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
    max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
    ``ADAPTIVE_POINT_CAPACITY_FACTOR``
    and `A` the capacity of the adaptive parameter search,
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`.
    As there, this is a **practical capacity, not a proven bound**: how many reference points a
    construction builds is only known once it has run, and a shuffled data set whose
    construction builds more fails with invalid input (it is then counted as a failed
    permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    :func:`tensor_omics.calc_js_comp_test_n_top_k_jsds` computes it.

    Rejects an out-of-range argument with invalid input (a non-finite significance level with
    `ERR_NAN_INF`), and, naming neither argument, the joint bounds on `max_n_genes_all_studies`
    and `n_studies` that `calc_adaptive_js_comp_test_bounds` rejects; the outputs are then left
    undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.
    k_start : int
        Pooled entries every neighborhood of the fixed-knob run takes unconditionally
        The minimum valid value is `1`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps of the parameter search, when it is rerun
        The minimum valid value is `1`.
        The default value is `1000`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Practical capacity for reference points
        n_bootstrapping_top_k_jsds : int
            Bootstrap heap size (values kept at each end of the bootstrap distribution)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_impl::calc_adaptive_js_comp_test_permutation_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    n_bootstrapping_top_k_jsds = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_adaptive_js_comp_test_permutation_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(k_start)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(n_bootstrapping_top_k_jsds),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_ADAPTIVE_JS_COMP_TEST_PERMUTATION_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "n_bootstrapping_top_k_jsds": n_bootstrapping_top_k_jsds.value,
    }
