#!/usr/bin/env python3
"""Detect subfunctionalization, neofunctionalization and dosage effects (SND).

After a duplication, the copies of an ancestral gene can go three ways, and
Tensor Omics tests for each separately:

  subfunctionalization  a SUBSET of copies sums back to the ancestor, each
                        copy having taken part of its expression
  dosage effect         copies point the same way as the ancestor but sum to
                        appreciably more of it
  neofunctionalization  a copy expresses in a direction the ancestor did not

The first two search over subsets of paralogs, so they work on raw expression
vectors, where magnitudes add up. The third asks only about direction, so it
works on RAP-projected unit vectors -- a different input space, and the most
common way to misuse this module.

Subsets are bitmasks: a result names its members as set bits, decoded with
mask_check_state. The ancestor in real use is the family's orthologs-only
centroid; see the gene family centroid recipe.

The family here is SYNTHETIC and built so each pattern is present exactly once,
which makes the output checkable. The repository ships no ortholog/paralog
designation, so no bundled family could serve.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/snd_detection.py
"""

import sys

import numpy as np

from tensor_omics import (calc_work_arr_paralog_subsets_size, detect_dosage_effect,
                          detect_neofunctionalization, detect_subfunctionalization,
                          filter_paralogs_by_pattern_dosage_effect,
                          filter_paralogs_by_pattern_subfunctionalization,
                          mask_check_state, mask_chunk_count, normalize_unit_length,
                          omics_vector_RAP_projection)


def angles_to(ancestor, genes):
    """Angle in radians between each gene and the ancestor, which the detectors need.

    No routine supplies these -- the filters take them as given, in [0, pi].
    """
    a_norm = np.linalg.norm(ancestor)
    cosines = [
        np.clip(gene @ ancestor / (np.linalg.norm(gene) * a_norm), -1.0, 1.0)
        for gene in genes.T
    ]
    return np.arccos(np.asarray(cosines))


def rap_unit(matrix):
    """RAP-project each column and scale it to unit length.

    What neofunctionalization expects, and what the other two must NOT be given.
    normalize_unit_length works in place and returns None, so the column is
    copied out, modified, and written back.
    """
    projected = np.asarray(omics_vector_RAP_projection(
        np.asfortranarray(matrix),
        np.ones(matrix.shape[1], dtype=bool),
        np.ones(matrix.shape[0], dtype=bool),
    ))
    out = np.empty(projected.shape, order="F")
    for i_col in range(projected.shape[1]):
        column = projected[:, i_col].copy()
        normalize_unit_length(column)
        out[:, i_col] = column
    return np.asfortranarray(out)


def members(mask_column, n_genes):
    """Decode one bitmask into the 1-based gene indices it contains."""
    return [i for i in range(1, n_genes + 1) if mask_check_state(mask_column, i)]


def make_family():
    """An ancestor and five paralogs, one pattern planted per group."""
    ancestor = np.array([6.0, 6.0, 6.0, 2.0])
    paralogs = {
        "sub_a": np.array([6.0, 6.0, 0.0, 1.0]),   # took two tissues
        "sub_b": np.array([0.0, 0.0, 6.0, 1.0]),   # took the others: sub_a+sub_b ~ ancestor
        "dos_a": np.array([5.0, 5.0, 5.0, 1.7]),   # same direction as the ancestor
        "dos_b": np.array([4.0, 4.0, 4.0, 1.3]),   # dos_a+dos_b clearly exceeds it
        "neo":   np.array([1.0, 1.0, 1.0, 9.0]),   # a domain the ancestor lacks
    }
    names = list(paralogs)
    return ancestor, names, np.asfortranarray(np.column_stack([paralogs[n] for n in names]))


def main(argv=None):
    ancestor, names, genes = make_family()
    n_genes = genes.shape[1]
    gene_to_fam = np.ones(n_genes, dtype=np.int32)
    n_chunks = mask_chunk_count(n_genes)

    angles = angles_to(ancestor, genes)
    print("family of %d paralogs, angle to ancestor (rad):" % n_genes)
    for name, angle in zip(names, angles):
        print(f"  {name:6} {angle:.4f}")
    median = float(np.median(angles))
    print(f"median angle {median:.4f}, used as the filter threshold\n")

    # --- subfunctionalization: wide-angle copies whose sum rebuilds the ancestor
    mask = np.asarray(filter_paralogs_by_pattern_subfunctionalization(
        angles, median, 1, gene_to_fam, n_chunks))[:, 0]
    print("subfunctionalization, kept by the filter:",
          [names[i - 1] for i in members(mask, n_genes)])
    sizing = calc_work_arr_paralog_subsets_size(n_genes, n_genes, mask)
    norms = np.linalg.norm(genes, axis=0)
    result = detect_subfunctionalization(
        ancestor, genes, mask, sizing["max_subset_size"],
        1.0,                                        # residual tolerated, the RDI threshold
        norms, (np.argsort(norms) + 1).astype(np.int32),
    )
    subsets = np.asarray(result["work_arr_paralog_subsets"])
    for k in range(result["n_results"]):
        picked = members(subsets[:, k], n_genes)
        residual = np.linalg.norm(ancestor - genes[:, [i - 1 for i in picked]].sum(axis=1))
        print(f"  subset {[names[i-1] for i in picked]}  residual {residual:.3f}")

    # --- dosage effect: narrow-angle copies summing to more than the ancestor
    mask = np.asarray(filter_paralogs_by_pattern_dosage_effect(
        angles, median, 1, gene_to_fam, n_chunks))[:, 0]
    print("\ndosage effect, kept by the filter:",
          [names[i - 1] for i in members(mask, n_genes)])
    sizing = calc_work_arr_paralog_subsets_size(n_genes, n_genes, mask)
    # max_angle defaults to pi, i.e. no angle constraint at all -- set it.
    result = detect_dosage_effect(ancestor, genes, mask, sizing["max_subset_size"],
                                  max_angle=0.3, gain_gamma=0.1)
    subsets = np.asarray(result["work_arr_paralog_subsets"])
    for k in range(result["n_results"]):
        picked = members(subsets[:, k], n_genes)
        gain = np.linalg.norm(genes[:, [i - 1 for i in picked]].sum(axis=1)) / np.linalg.norm(ancestor)
        print(f"  subset {[names[i-1] for i in picked]}  magnitude {gain:.2f}x the ancestor")

    # --- neofunctionalization: direction only, so RAP unit vectors
    ancestors_u = rap_unit(ancestor.reshape(-1, 1))
    genes_u = rap_unit(genes)
    neofunc = np.asarray(detect_neofunctionalization(
        ancestors_u, genes_u, gene_to_fam, np.full(genes.shape[0], 0.5)))
    print("\nneofunctionalization, axes gained or lost beyond 0.5:")
    for i_gene, name in enumerate(names):
        axes = [a + 1 for a in range(neofunc.shape[1]) if neofunc[i_gene, a]]
        print(f"  {name:6} {'axes ' + str(axes) if axes else '-'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
