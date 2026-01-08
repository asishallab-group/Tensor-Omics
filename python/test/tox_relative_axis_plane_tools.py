import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions import tox_vector_RAP_projection, tox_field_RAP_projection


def test_omics_vector_RAP_projection_call():
    """
    Test the tox_vector_RAP_projection function for correct projection behavior.

    Generates random vectors and selection masks, computes projections, and asserts that the sum of each projected column is zero.
    Prints a success message if the test passes.
    """
    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    vecs = np.random.rand(n_axes, n_vecs)
    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    projections = tox_vector_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

    for i_vec in range(projections.shape[1]):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"
    print("test_omics_vector_RAP_projection_call PASSED")


def test_omics_field_RAP_projection_call():
    """
    Test the tox_field_RAP_projection function for correct field projection behavior.

    Generates random field vectors and selection masks, computes projections, and asserts that the sum of each projected column is zero.
    Prints a success message if the test passes.
    """
    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    vecs = np.random.rand(2 * n_axes, n_vecs)
    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    projections = tox_field_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

    for i_vec in range(projections.shape[1]):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"
    print("test_omics_field_RAP_projection_call PASSED")


def main():
    """
    Main function to run all RAP projection tests.

    Runs both vector and field RAP projection tests and prints a summary message if all pass.
    """
    test_omics_vector_RAP_projection_call()
    test_omics_field_RAP_projection_call()
    print("All RAP projection tests PASSED.")


if __name__ == '__main__':
    """
    Entry point for running the RAP projection test suite.

    Calls the main() function to execute all tests.
    """
    main()
