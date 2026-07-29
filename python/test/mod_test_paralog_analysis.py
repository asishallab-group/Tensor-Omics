import numpy as np
import math
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    mask_check_state,
    filter_paralogs_by_pattern_dosage_effect,
    filter_paralogs_by_pattern_subfunctionalization,
    calc_work_arr_paralog_subsets_size,
    detect_dosage_effect,
    detect_subfunctionalization,
    mask_chunk_count,
    detect_neofunctionalization,
    normalize_unit_length
)
from test_helpers import run_all_tests


def _normalize_columns(matrix):
    """Column-wise unit normalization; the Fortran routine works in place on a vector."""
    out = np.array(matrix, dtype=np.float64, order="F")
    for i_col in range(out.shape[1]):
        column = np.ascontiguousarray(out[:, i_col])
        normalize_unit_length(column)
        out[:, i_col] = column
    return out


def test_paralog_functions():
    # Testing Mask Logic

    n_paralogs = 5
    i_paralog = 2
    chunk_count = mask_chunk_count(n_paralogs)
    assert chunk_count == 1, f"Chunk count for {n_paralogs} paralogs should be 1, got " + str(chunk_count)

    bit_mask = np.zeros(chunk_count, dtype=np.int32)
    bit_mask[i_paralog // 32] = 1 << (i_paralog % 32)
    state = mask_check_state(bit_mask, i_paralog + 1)
    assert state, f"Paralog {i_paralog + 1} should be active"

    # Testing Pattern Filtering

    n_families = 1
    gene_to_fam = np.full(n_paralogs, 1)
    angles = np.array([0.1, 0.3, 0.5, 0.7, 0.9], dtype=np.float64)
    threshold = 0.6

    gene_to_fam = gene_to_fam.astype(np.int32)
    # masks come back as (n_mask_chunks, n_families); the detectors take one family's column
    dosage_mask = filter_paralogs_by_pattern_dosage_effect(
        angles, threshold, n_families, gene_to_fam, chunk_count)[:, 0]
    assert (dosage_mask == [7]).all(), "dosage_mask should be [7], got " + str(dosage_mask)

    subfunc_mask = filter_paralogs_by_pattern_subfunctionalization(
        angles, threshold, n_families, gene_to_fam, chunk_count)[:, 0]
    # 24 == 0b11000: genes 4 and 5, the ones whose angle exceeds the threshold
    assert (subfunc_mask == [24]).all(), "subfunc_mask should be [24], got " + str(subfunc_mask)

    # Testing Work Array Size Calculation

    max_subset_size = 3
    sizes = calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, dosage_mask)
    work_array_size = sizes["work_array_size"]
    actual_max_subset_size = sizes["max_subset_size"]
    assert work_array_size == 3, "Expected work_array_size to be 3, got " + str(work_array_size)
    assert actual_max_subset_size == 3, "Expected actual_max_subset_size to be 3, got " + str(actual_max_subset_size)

    # Testing Dosage Effect Detection

    ancestor = np.array([1.0, 1.0], dtype=np.float64, order="F")
    paralogs = np.array([
        [1.1, 1.2, 1.3, 1.4, 1.5],
        [0.9, 0.8, 0.7, 0.6, 0.5]
    ], dtype=np.float64, order="F")

    dosage_result = detect_dosage_effect(
        ancestor=ancestor,
        genes=paralogs,
        filtered_paralogs_mask=dosage_mask,
        max_subset_size=actual_max_subset_size,
        gain_gamma=0.1,
        max_angle=math.pi
    )
    assert (dosage_result['work_arr_paralog_subsets'] == [[6, 3, 5]]).all(), "Dosage Effect results should be [[6, 3, 5]], got " + str(dosage_result['work_arr_paralog_subsets'])
    assert dosage_result['n_results'] == 3, "Expected Dosage Effect n_results to be 3, got " + str(dosage_result['n_results'])

    # an over-large max_subset_size is capped by calc_work_arr_paralog_subsets_size (inout
    # feedback): the same result comes back at once, not after a runaway subset-size loop
    capped = detect_dosage_effect(
        ancestor=ancestor,
        genes=paralogs,
        filtered_paralogs_mask=dosage_mask,
        max_subset_size=2_000_000_000,
        gain_gamma=0.1,
        max_angle=math.pi
    )
    assert capped['n_results'] == 3, "Over-large max_subset_size should be capped, got n_results " + str(capped['n_results'])
    assert (capped['work_arr_paralog_subsets'] == [[6, 3, 5]]).all(), "Capped run should match, got " + str(capped['work_arr_paralog_subsets'])

    # Testing Subfunctionalization Detection

    norms = np.sqrt(np.sum(paralogs**2, axis=0))
    sorted_perm = np.argsort(norms).astype(np.int32) + 1

    subfunc_sizes = calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, subfunc_mask)
    subfunc_work_size = subfunc_sizes["work_array_size"]
    subfunc_max_subset_size = subfunc_sizes["max_subset_size"]

    subfunc_result = detect_subfunctionalization(
        ancestor=ancestor,
        genes=paralogs,
        rdi_threshold=0.5,
        filtered_paralogs_mask=subfunc_mask,
        max_subset_size=subfunc_max_subset_size,
        paralog_norms=norms,
        sorted_paralog_norms_perm=sorted_perm
    )
    assert (subfunc_result['work_arr_paralog_subsets'] == []).all(), "Subfunctionalization results should be empty array, got " + str(subfunc_result['work_arr_paralog_subsets'])
    assert subfunc_result['n_results'] == 0, "Expected Subfunctionalization n_results to be zero, got " + str(subfunc_result['n_results'])

    # Testing Edge Cases

    try:
        empty_mask = mask_chunk_count(0)
    except Exception as e:
        raise AssertionError("mask_chunk_count throws error for empty mask")

    single_mask = mask_chunk_count(1)
    assert single_mask == 1, "Single paralog mask chunk count should be 1, got " + str(single_mask)


def test_detect_neofunctionalization():
    # -------------------------------
    # Case 1: Differences below threshold → all false (all zeros)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # normalize_unit_length works in place on a contiguous vector
    ancestors = _normalize_columns(ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.05, 0.05], dtype=np.float64, order="F")

    # Build genes identical to ancestors for each gene's family
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1]  # adjust index for Python 0-based

    neofunc = detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.zeros((3, 2), dtype=np.int32, order="F")
    assert np.array_equal(neofunc, expected), "Case 1 output mismatch"

    # -------------------------------
    # Case 2: Differences above threshold → some true (some ones)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # normalize_unit_length works in place on a contiguous vector
    ancestors = _normalize_columns(ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.2, 0.2], dtype=np.float64, order="F")

    # Build genes offset by threshold * family index
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1] - thresholds * gene_to_fam[i_gene]

    neofunc = detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.array([[False, False],
                         [True, True],
                         [False, False]], dtype=bool, order="F")
    assert np.array_equal(neofunc, expected), "Case 2 output mismatch"


if __name__ == "__main__":
    run_all_tests(globals().values())
