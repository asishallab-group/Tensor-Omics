import numpy as np
import math
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_mask_check_state,
    tox_filter_paralogs_by_pattern_dosage_effect,
    tox_filter_paralogs_by_pattern_subfunctionalization,
    tox_calc_work_arr_paralog_subsets_size,
    tox_detect_dosage_effect,
    tox_detect_subfunctionalization,
    tox_mask_chunk_count,
    tox_detect_neofunctionalization,
    tox_normalize_unit_length
)


def test_paralog_functions():
    # Testing Mask Logic

    n_paralogs = 5
    i_paralog = 2
    chunk_count = tox_mask_chunk_count(n_paralogs)
    assert chunk_count == 1, f"Chunk count for {n_paralogs} paralogs should be 1, got " + str(chunk_count)

    bit_mask = np.zeros(chunk_count, dtype=np.int32)
    bit_mask[i_paralog // 32] = 1 << (i_paralog % 32)
    state = tox_mask_check_state(bit_mask, i_paralog + 1)
    assert state, f"Paralog {i_paralog + 1} should be active"

    # Testing Pattern Filtering

    n_families = 1
    gene_to_fam = np.full(n_paralogs, 1)
    angles = np.array([0.1, 0.3, 0.5, 0.7, 0.9], dtype=np.float64)
    threshold = 0.6

    dosage_mask, = tox_filter_paralogs_by_pattern_dosage_effect(angles, threshold, gene_to_fam, n_families)
    assert (dosage_mask == [7]).all(), "dosage_mask should be [7], got " + str(dosage_mask)

    subfunc_mask, = tox_filter_paralogs_by_pattern_subfunctionalization(angles, threshold, gene_to_fam, n_families)
    assert (subfunc_mask == []).all(), "subfunc_mask should be [], got " + str(subfunc_mask)

    # Testing Work Array Size Calculation

    max_subset_size = 3
    work_array_size, actual_max_subset_size = tox_calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, dosage_mask).values()
    assert work_array_size == 3, "Expected work_array_size to be 3, got " + str(work_array_size)
    assert actual_max_subset_size == 3, "Expected actual_max_subset_size to be 3, got " + str(actual_max_subset_size)

    # Testing Dosage Effect Detection

    ancestor = np.array([1.0, 1.0], dtype=np.float64, order="F")
    paralogs = np.array([
        [1.1, 1.2, 1.3, 1.4, 1.5],
        [0.9, 0.8, 0.7, 0.6, 0.5]
    ], dtype=np.float64, order="F")

    dosage_result = tox_detect_dosage_effect(
        ancestor=ancestor,
        genes=paralogs,
        filtered_paralogs_mask=dosage_mask,
        max_subset_size=n_paralogs,
        gain_gamma=0.1,
        max_angle=math.pi
    )
    assert (dosage_result['results'] == [[6, 3, 5]]).all(), "Dosage Effect results should be [[6, 3, 5]], got " + str(dosage_result['results'])
    assert dosage_result['n_results'] == 3, "Expected Dosage Effect n_results to be 3, got " + str(dosage_result['n_results'])

    # Testing Subfunctionalization Detection

    norms = np.sqrt(np.sum(paralogs**2, axis=0))
    sorted_perm = np.argsort(norms).astype(np.int32) + 1

    subfunc_result = tox_detect_subfunctionalization(
        ancestor=ancestor,
        genes=paralogs,
        rdi_threshold=0.5,
        filtered_paralogs_mask=subfunc_mask,
        max_subset_size=n_paralogs,
        paralog_norms=norms,
        sorted_paralog_norms_perm=sorted_perm
    )
    assert (subfunc_result['results'] == []).all(), "Subfunctionalization results should be empty array, got " + str(subfunc_result['results'])
    assert subfunc_result['n_results'] == 0, "Expected Subfunctionalization n_results to be zero, got " + str(subfunc_result['n_results'])

    # Testing Edge Cases

    try:
        empty_mask = tox_mask_chunk_count(0)
    except Exception as e:
        raise AssertionError("tox_mask_chunk_count throws error for empty mask")

    single_mask = tox_mask_chunk_count(1)
    assert single_mask == 1, "Single paralog mask chunk count should be 1, got " + str(single_mask)

    print("✅ Paralog functions passed.")


def test_detect_neofunctionalization():
    # -------------------------------
    # Case 1: Differences below threshold → all false (all zeros)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # Normalize ancestors
    ancestors = np.apply_along_axis(tox_normalize_unit_length, axis=0, arr=ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.05, 0.05], dtype=np.float64, order="F")

    # Build genes identical to ancestors for each gene's family
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1]  # adjust index for Python 0-based

    neofunc = tox_detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.zeros((3, 2), dtype=np.int32, order="F")
    assert np.array_equal(neofunc, expected), "Case 1 output mismatch"

    # -------------------------------
    # Case 2: Differences above threshold → some true (some ones)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # Normalize ancestors
    ancestors = np.apply_along_axis(tox_normalize_unit_length, axis=0, arr=ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.2, 0.2], dtype=np.float64, order="F")

    # Build genes offset by threshold * family index
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1] - thresholds * gene_to_fam[i_gene]

    neofunc = tox_detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.array([[False, False],
                         [True, True],
                         [False, False]], dtype=bool, order="F")
    assert np.array_equal(neofunc, expected), "Case 2 output mismatch"

    print("✅ Neofunctionalization passed.")


def main():
    print("=================================================")
    print("    TOX PARALOG ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    # Run the tests
    test_paralog_functions()
    test_detect_neofunctionalization()


if __name__ == '__main__':
    main()
