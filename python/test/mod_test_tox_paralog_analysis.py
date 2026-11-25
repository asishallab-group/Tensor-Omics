"""
The `tox_paralog_analysis` suite: `detect_neofunctionalization`, `detect_dosage_effect`,
`detect_subfunctionalization`, the two `filter_paralogs_by_pattern_*` filters and
`calc_paralog_subsets_bitmasks_size` can be called from Python, return the documented types and
shapes, and raise the documented errors, naming the argument. Whether the subsets and flags are
right is the Fortran suite's job (test/mod_test_tox_paralog_analysis/), as the coding guide's
"Where to Test What" asks. The few values checked here are the binding's own: genes are bits in
index order, from 1, a bit mask word with bit 31 set comes back as the int32 -2**31, and the
binding sizes the work array and cuts max_subset_size itself.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (detect_neofunctionalization, detect_dosage_effect, detect_subfunctionalization,
                          filter_paralogs_by_pattern_dosage_effect, filter_paralogs_by_pattern_subfunctionalization,
                          calc_paralog_subsets_bitmasks_size)
from test_helpers import run_all_tests, assert_error, assert_refused_by_name
from tensor_omics.error_handling import ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF

#: gene_to_fam of a gene in no family (detect_neofunctionalization only)
UNASSIGNED = 0
#: the int32 word whose only set bit is bit 31, gene 32 of a mask
ONLY_BIT_31 = -2**31


def _neofunc_inputs():
    """Ancestors (axes, families) [0.5, 0.5] and [-0.5, 0]; three genes, in families 2, 1 and none;
    thresholds 0.25."""
    ancestors = np.array([[0.5, 0.5], [-0.5, 0.0]]).T.copy(order='F')
    genes = np.array([[0.5, 0.5], [0.5, 0.0], [-1.0, 1.0]]).T.copy(order='F')
    gene_to_fam = np.array([2, 1, UNASSIGNED], dtype=np.int32)
    thresholds = np.array([0.25, 0.25])
    return ancestors, genes, gene_to_fam, thresholds


def _pattern_inputs():
    """Three genes, each 0.75 times the ancestor [1, 1], all kept (the word 7): every pair is 1.5
    times the ancestor, a dosage effect."""
    ancestor = np.array([1.0, 1.0])
    genes = np.full((2, 3), 0.75, order='F')
    return ancestor, genes, np.array([7], dtype=np.int32)


def _subfunc_inputs():
    """Ancestor [2, 0], three genes [1, 0], all kept, with their norms and ascending order."""
    ancestor = np.array([2.0, 0.0])
    genes = np.array([[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]], order='F')
    return ancestor, genes, np.array([7], dtype=np.int32), np.ones(3), np.array([1, 2, 3], dtype=np.int32)


def _filter_inputs():
    """Angles [0, 1, 0], one family."""
    return np.array([0.0, 1.0, 0.0]), np.ones(3, dtype=np.int32)


# =====================
# detect_neofunctionalization
# =====================

def test_detect_neofunctionalization():
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    result = detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    assert isinstance(result, np.ndarray), f"expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.bool_, f"expected bool, got {result.dtype}"
    assert result.shape == (3, 2), f"expected (genes, axes) = (3, 2), got {result.shape}"
    assert not result.flags.writeable, "a result is a value, it should be read-only"


def test_detect_neofunctionalization_takes_family_ids_1_based():
    # the one value: id 2 is the second column of ancestors, and 0 has none. Gene 1 is family 1's
    # ancestor exactly, so it is flagged only if compared with the second.
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    result = detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.array([[True, True], [False, True], [False, False]])
    assert np.array_equal(result, expected), f"got {result.tolist()}"


def test_detect_neofunctionalization_rejects_a_family_id_past_n_families():
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    gene_to_fam[-1] = 3
    assert_error(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds),
                 "a family id past n_families", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_detect_neofunctionalization_rejects_genes_out_of_range():
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    genes[-1, -1] = 2.0
    assert_error(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds),
                 "a gene component above 1", ERR_INVALID_INPUT, argument="genes")


def test_detect_neofunctionalization_rejects_nan_ancestor():
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    ancestors[-1, -1] = np.nan
    assert_error(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds),
                 "a NaN in ancestors", ERR_NAN_INF, argument="ancestors")


def test_detect_neofunctionalization_rejects_infinite_threshold():
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    thresholds[-1] = np.inf
    assert_error(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds),
                 "an infinite threshold", ERR_NAN_INF, argument="thresholds")


def test_detect_neofunctionalization_rejects_no_families():
    # n_families is read off ancestors
    _, genes, _, thresholds = _neofunc_inputs()
    assert_error(lambda: detect_neofunctionalization(np.empty((2, 0), order='F'), genes,
                                                     np.zeros(3, dtype=np.int32), thresholds),
                 "no families", ERR_EMPTY_INPUT, argument="ancestors")


def test_detect_neofunctionalization_rejects_no_genes():
    # n_genes is read off genes
    ancestors, _, _, thresholds = _neofunc_inputs()
    assert_error(lambda: detect_neofunctionalization(ancestors, np.empty((2, 0), order='F'),
                                                     np.empty(0, dtype=np.int32), thresholds),
                 "no genes", ERR_EMPTY_INPUT, argument="genes")


def test_detect_neofunctionalization_rejects_mismatched_shapes():
    # the binding's own shape checks, before the library is called
    ancestors, genes, gene_to_fam, thresholds = _neofunc_inputs()
    assert_refused_by_name(lambda: detect_neofunctionalization(ancestors, np.zeros((3, 3), order='F'), gene_to_fam,
                                                               thresholds),
                           "genes over another number of axes", "genes")
    assert_refused_by_name(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam[:2], thresholds),
                           "a mapping of another number of genes", "gene_to_fam")
    assert_refused_by_name(lambda: detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds[:1]),
                           "thresholds of another number of axes", "thresholds")


# =====================
# detect_dosage_effect
# =====================

def test_detect_dosage_effect():
    ancestor, genes, mask = _pattern_inputs()
    result = detect_dosage_effect(ancestor, genes, mask, 3)
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == {"n_results", "paralog_subsets_bitmasks"}, f"got the keys {sorted(result)}"
    assert isinstance(result["n_results"], int), f"expected an int, got {type(result['n_results']).__name__}"
    subsets = result["paralog_subsets_bitmasks"]
    assert subsets.dtype == np.int32, f"expected int32 words, got {subsets.dtype}"
    # the whole work array comes back, as many subsets as the sizing asks for
    work_array_size = calc_paralog_subsets_bitmasks_size(3, 3, mask)["work_array_size"]
    assert subsets.shape == (1, work_array_size), f"expected (n_mask_chunks, {work_array_size}), got {subsets.shape}"
    assert not subsets.flags.writeable, "a result is a value, it should be read-only"


def test_detect_dosage_effect_numbers_genes_from_bit_0():
    # the one value: the three pairs are genes {1, 2}, {1, 3} and {2, 3}, the words 3, 5 and 6
    ancestor, genes, mask = _pattern_inputs()
    result = detect_dosage_effect(ancestor, genes, mask, 3)
    found = sorted(result["paralog_subsets_bitmasks"][0, :result["n_results"]].tolist())
    assert found == [3, 5, 6], f"got {found}"


def test_detect_dosage_effect_cuts_max_subset_size():
    # the binding sizes the work array itself, and the sizing cuts an over-large max_subset_size
    # to the kept genes, so a search size of 2e9 asks for no more than three genes do
    ancestor, genes, mask = _pattern_inputs()
    result = detect_dosage_effect(ancestor, genes, mask, 2_000_000_000)
    assert result["paralog_subsets_bitmasks"].shape == (1, 3), f"got {result['paralog_subsets_bitmasks'].shape}"
    assert result["n_results"] == 3, f"got {result['n_results']}"


def test_detect_dosage_effect_returns_bit_31_as_a_negative_word():
    # genes 32 and 33 are bit 31 of the first word and bit 0 of the second; the pair comes back
    # as the same two words
    ancestor = np.array([1.0, 1.0])
    genes = np.full((2, 33), 0.75, order='F')
    mask = np.array([ONLY_BIT_31, 1], dtype=np.int32)
    result = detect_dosage_effect(ancestor, genes, mask, 2)
    assert result["n_results"] == 1, f"got {result['n_results']}"
    assert result["paralog_subsets_bitmasks"][:, 0].tolist() == [ONLY_BIT_31, 1], \
        f"got {result['paralog_subsets_bitmasks'][:, 0].tolist()}"


def test_detect_dosage_effect_with_fewer_than_two_kept_genes():
    # BUG (fails today): a family whose filter kept no gene, or one, holds no pair, and the
    # detection has nothing to find. The binding sizes the work array at 0 subsets and passes
    # that on, which the library refuses as ERR_EMPTY_INPUT naming n_paralog_subsets -- an
    # argument the caller cannot pass. Expected: no result, n_results 0.
    ancestor, genes, _ = _pattern_inputs()
    for kept in (0, 2):
        result = detect_dosage_effect(ancestor, genes, np.array([kept], dtype=np.int32), 3)
        assert result["n_results"] == 0, f"the word {kept}: got {result['n_results']}"


def test_detect_dosage_effect_rejects_a_negative_max_subset_size():
    ancestor, genes, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(ancestor, genes, mask, -1),
                 "max_subset_size -1", ERR_INVALID_INPUT, argument="max_subset_size")


def test_detect_dosage_effect_rejects_a_mask_too_short_for_the_genes():
    # 40 genes need two words; the sizing runs first and blames the genes it was handed
    ancestor, _, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(ancestor, np.ones((2, 40), order='F'), mask, 3),
                 "one word for 40 genes", ERR_INVALID_INPUT, argument="genes")


def test_detect_dosage_effect_rejects_a_gain_of_zero():
    ancestor, genes, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(ancestor, genes, mask, 3, gain_gamma=0.0),
                 "gain_gamma 0", ERR_INVALID_INPUT, argument="gain_gamma")


def test_detect_dosage_effect_rejects_an_angle_above_pi():
    ancestor, genes, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(ancestor, genes, mask, 3, max_angle=4.0),
                 "max_angle 4", ERR_INVALID_INPUT, argument="max_angle")


def test_detect_dosage_effect_rejects_nan_ancestor():
    _, genes, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(np.array([1.0, np.nan]), genes, mask, 3),
                 "a NaN in ancestor", ERR_NAN_INF, argument="ancestor")


def test_detect_dosage_effect_rejects_infinite_genes():
    ancestor, genes, mask = _pattern_inputs()
    genes[-1, -1] = -np.inf
    assert_error(lambda: detect_dosage_effect(ancestor, genes, mask, 3),
                 "an infinite gene", ERR_NAN_INF, argument="genes")


def test_detect_dosage_effect_rejects_a_zero_ancestor():
    _, genes, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(np.zeros(2), genes, mask, 3),
                 "an ancestor of zeros", ERR_INVALID_INPUT, argument="ancestor")


def test_detect_dosage_effect_rejects_no_dimensions():
    # n_dims is read off ancestor
    _, _, mask = _pattern_inputs()
    assert_error(lambda: detect_dosage_effect(np.array([]), np.empty((0, 3), order='F'), mask, 3),
                 "no dimensions", ERR_EMPTY_INPUT, argument="ancestor")


def test_detect_dosage_effect_rejects_genes_of_other_dimensions():
    # the binding's own shape check, before the library is called
    _, genes, mask = _pattern_inputs()
    assert_refused_by_name(lambda: detect_dosage_effect(np.ones(3), genes, mask, 3),
                           "genes over another number of dimensions", "genes")


# =====================
# detect_subfunctionalization
# =====================

def test_detect_subfunctionalization():
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    result = detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending)
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == {"n_results", "paralog_subsets_bitmasks"}, f"got the keys {sorted(result)}"
    assert isinstance(result["n_results"], int), f"expected an int, got {type(result['n_results']).__name__}"
    subsets = result["paralog_subsets_bitmasks"]
    assert subsets.dtype == np.int32, f"expected int32 words, got {subsets.dtype}"
    work_array_size = calc_paralog_subsets_bitmasks_size(3, 3, mask)["work_array_size"]
    assert subsets.shape == (1, work_array_size), f"expected (n_mask_chunks, {work_array_size}), got {subsets.shape}"
    assert not subsets.flags.writeable, "a result is a value, it should be read-only"


def test_detect_subfunctionalization_numbers_genes_from_bit_0():
    # the one value: any pair of the three [1, 0] genes rebuilds [2, 0], the words 3, 5 and 6
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    result = detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending)
    found = sorted(result["paralog_subsets_bitmasks"][0, :result["n_results"]].tolist())
    assert found == [3, 5, 6], f"got {found}"


def test_detect_subfunctionalization_with_fewer_than_two_kept_genes():
    # BUG (fails today): as test_detect_dosage_effect_with_fewer_than_two_kept_genes, refused as
    # ERR_EMPTY_INPUT naming n_paralog_subsets. Expected: no result, n_results 0.
    ancestor, genes, _, norms, ascending = _subfunc_inputs()
    for kept in (0, 2):
        result = detect_subfunctionalization(ancestor, genes, np.array([kept], dtype=np.int32), 3, 0.0, norms, ascending)
        assert result["n_results"] == 0, f"the word {kept}: got {result['n_results']}"


def test_detect_subfunctionalization_rejects_a_negative_threshold():
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    assert_error(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, -1.0, norms, ascending),
                 "rdi_threshold -1", ERR_INVALID_INPUT, argument="rdi_threshold")


def test_detect_subfunctionalization_rejects_a_negative_norm():
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    norms[-1] = -1.0
    assert_error(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending),
                 "a norm of -1", ERR_INVALID_INPUT, argument="paralog_norms")


def test_detect_subfunctionalization_rejects_a_permutation_entry_of_0():
    # the entries are 1-based gene indices
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    ascending[0] = 0
    assert_error(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending),
                 "a permutation entry 0", ERR_INVALID_INPUT, argument="sorted_paralog_norms_perm")


def test_detect_subfunctionalization_rejects_a_nan_norm():
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    norms[-1] = np.nan
    assert_error(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending),
                 "a NaN norm", ERR_NAN_INF, argument="paralog_norms")


def test_detect_subfunctionalization_rejects_mismatched_shapes():
    # the binding's own shape checks, before the library is called
    ancestor, genes, mask, norms, ascending = _subfunc_inputs()
    assert_refused_by_name(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms[:2], ascending),
                           "norms of another number of genes", "paralog_norms")
    assert_refused_by_name(lambda: detect_subfunctionalization(ancestor, genes, mask, 3, 0.0, norms, ascending[:2]),
                           "a permutation of another number of genes", "sorted_paralog_norms_perm")
    assert_refused_by_name(lambda: detect_subfunctionalization(np.ones(3), genes, mask, 3, 0.0, norms, ascending),
                           "genes over another number of dimensions", "genes")


# =====================
# filter_paralogs_by_pattern_dosage_effect / _subfunctionalization
# =====================

FILTERS = (filter_paralogs_by_pattern_dosage_effect, filter_paralogs_by_pattern_subfunctionalization)


def test_filter_paralogs_by_pattern():
    angles, gene_to_fam = _filter_inputs()
    for filter_paralogs in FILTERS:
        result = filter_paralogs(angles, 0.5, 1, gene_to_fam, 1)
        assert isinstance(result, np.ndarray), f"{filter_paralogs.__name__}: got {type(result).__name__}"
        assert result.dtype == np.int32, f"{filter_paralogs.__name__}: expected int32 words, got {result.dtype}"
        assert result.shape == (1, 1), f"{filter_paralogs.__name__}: expected (n_mask_chunks, n_families), got {result.shape}"
        assert not result.flags.writeable, f"{filter_paralogs.__name__}: a result is a value, it should be read-only"


def test_filter_paralogs_by_pattern_numbers_genes_from_bit_0():
    # the one value: angles [0, 1, 0] at 0.5 keep genes 1 and 3 for dosage effect (the word 5),
    # gene 2 for subfunctionalization (the word 2)
    angles, gene_to_fam = _filter_inputs()
    dosage = filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1, gene_to_fam, 1)
    subfunc = filter_paralogs_by_pattern_subfunctionalization(angles, 0.5, 1, gene_to_fam, 1)
    assert dosage.tolist() == [[5]] and subfunc.tolist() == [[2]], f"got {dosage.tolist()} and {subfunc.tolist()}"


def test_filter_paralogs_by_pattern_returns_bit_31_as_a_negative_word():
    # gene 32 alone is bit 31 of the word
    angles = np.ones(32)
    angles[-1] = 0.0
    result = filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1, np.ones(32, dtype=np.int32), 1)
    assert result.tolist() == [[ONLY_BIT_31]], f"got {result.tolist()}"


def test_filter_paralogs_by_pattern_rejects_more_families_than_genes():
    angles, gene_to_fam = _filter_inputs()
    for filter_paralogs in FILTERS:
        assert_error(lambda: filter_paralogs(angles, 0.5, 4, gene_to_fam, 1),
                     f"{filter_paralogs.__name__}: 4 families of 3 genes", ERR_INVALID_INPUT, argument="n_families")


def test_filter_paralogs_by_pattern_rejects_negative_extents():
    # the binding's own refusal, before it allocates the output by them
    angles, gene_to_fam = _filter_inputs()
    for filter_paralogs in FILTERS:
        assert_refused_by_name(lambda: filter_paralogs(angles, 0.5, -1, gene_to_fam, 1),
                               f"{filter_paralogs.__name__}: n_families -1", "n_families")
        assert_refused_by_name(lambda: filter_paralogs(angles, 0.5, 1, gene_to_fam, -1),
                               f"{filter_paralogs.__name__}: n_mask_chunks -1", "n_mask_chunks")


def test_filter_paralogs_by_pattern_rejects_too_few_words():
    angles, gene_to_fam = _filter_inputs()
    assert_error(lambda: filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1, gene_to_fam, 0),
                 "no word for 3 genes", ERR_INVALID_INPUT, argument="n_mask_chunks")


def test_filter_paralogs_by_pattern_rejects_a_family_id_of_0():
    # these filters take no unassigned genes
    angles, gene_to_fam = _filter_inputs()
    gene_to_fam[1] = 0
    assert_error(lambda: filter_paralogs_by_pattern_subfunctionalization(angles, 0.5, 1, gene_to_fam, 1),
                 "a family id 0", ERR_INVALID_INPUT, argument="gene_to_fam")


def test_filter_paralogs_by_pattern_rejects_an_angle_above_pi():
    angles, gene_to_fam = _filter_inputs()
    angles[1] = 4.0
    assert_error(lambda: filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1, gene_to_fam, 1),
                 "an angle of 4", ERR_INVALID_INPUT, argument="gene_angles")


def test_filter_paralogs_by_pattern_rejects_a_nan_threshold():
    angles, gene_to_fam = _filter_inputs()
    assert_error(lambda: filter_paralogs_by_pattern_subfunctionalization(angles, np.nan, 1, gene_to_fam, 1),
                 "a NaN threshold", ERR_NAN_INF, argument="threshold")


def test_filter_paralogs_by_pattern_rejects_no_genes():
    # n_genes is read off gene_angles
    assert_error(lambda: filter_paralogs_by_pattern_dosage_effect(np.array([]), 0.5, 1, np.empty(0, dtype=np.int32), 1),
                 "no genes", ERR_EMPTY_INPUT, argument="gene_angles")


def test_filter_paralogs_by_pattern_rejects_a_mapping_of_other_genes():
    # the binding's own shape check, before the library is called
    angles, gene_to_fam = _filter_inputs()
    assert_refused_by_name(lambda: filter_paralogs_by_pattern_dosage_effect(angles, 0.5, 1, gene_to_fam[:2], 1),
                           "a mapping of another number of genes", "gene_to_fam")


# =====================
# calc_paralog_subsets_bitmasks_size
# =====================

def test_calc_paralog_subsets_bitmasks_size():
    result = calc_paralog_subsets_bitmasks_size(3, 3, np.array([7], dtype=np.int32))
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == {"max_subset_size", "work_array_size"}, f"got the keys {sorted(result)}"
    assert all(isinstance(value, int) for value in result.values()), f"expected ints, got {result}"


def test_calc_paralog_subsets_bitmasks_size_hands_back_the_cut_size():
    # the one value: max_subset_size is in and out; 100 comes back as the three kept genes
    result = calc_paralog_subsets_bitmasks_size(100, 3, np.array([7], dtype=np.int32))
    assert result["max_subset_size"] == 3, f"got {result}"


def test_calc_paralog_subsets_bitmasks_size_rejects_a_negative_max_subset_size():
    assert_error(lambda: calc_paralog_subsets_bitmasks_size(-1, 3, np.array([7], dtype=np.int32)),
                 "max_subset_size -1", ERR_INVALID_INPUT, argument="max_subset_size")


def test_calc_paralog_subsets_bitmasks_size_rejects_more_genes_than_bits():
    assert_error(lambda: calc_paralog_subsets_bitmasks_size(3, 33, np.array([7], dtype=np.int32)),
                 "33 genes in one word", ERR_INVALID_INPUT, argument="n_genes")


def test_calc_paralog_subsets_bitmasks_size_rejects_no_genes():
    assert_error(lambda: calc_paralog_subsets_bitmasks_size(3, 0, np.array([7], dtype=np.int32)),
                 "n_genes 0", ERR_INVALID_INPUT, argument="n_genes")


if __name__ == "__main__":
    run_all_tests(globals().values())
