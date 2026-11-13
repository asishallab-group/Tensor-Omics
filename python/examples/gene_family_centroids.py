#!/usr/bin/env python3
"""Compute one expression centroid per gene family.

A family centroid is the element-wise mean of its members' expression vectors.
It is the reference every later Tensor Omics step measures against: distances
to it, shift vectors from it, and the outlier test built on those.

Two centroids are available and they answer different questions:

    group_centroid_all         mean over every gene in the family
    group_centroid_orthologs   mean over the family's orthologs only

Pass --orthologs to use the second. The ortholog assignment is *your* data --
it comes from an orthology pipeline (OrthoFinder, synteny/Cactus, OMAmer),
not from Tensor Omics -- so this script reads it as a plain list of gene ids.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/gene_family_centroids.py

which computes all-genes centroids for the bundled 458 families over
``material/normalization.tsv`` and writes ``results/family_centroids.tsv``.

Reading the input files is plumbing and lives in ``example_data.py``; what is
below is the workflow itself.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from example_data import (UNASSIGNED, build_gene_to_family, read_expression,
                          read_families, read_ortholog_mask, write_labeled_matrix)
from tensor_omics import group_centroid_all, group_centroid_orthologs

DEFAULT_FAMILIES = Path("material/filtered_families.tsv")
DEFAULT_EXPRESSION = Path("material/normalization.tsv")
DEFAULT_OUTPUT = Path("results/family_centroids.tsv")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES,
                        help=f"family table (default: {DEFAULT_FAMILIES})")
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION,
                        help=f"normalized gene x axis table (default: {DEFAULT_EXPRESSION})")
    parser.add_argument("--orthologs", type=Path, default=None,
                        help="gene ids of the orthologs, one per line; if given, "
                             "centroids are taken over these genes only")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT,
                        help=f"where to write the centroids (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args(argv)

    family_ids, genes_of_family = read_families(args.families)
    gene_ids, axes, expression_vectors = read_expression(args.expression)
    gene_to_family, n_missing = build_gene_to_family(gene_ids, genes_of_family)

    n_families = len(family_ids)
    n_assigned = int(np.count_nonzero(gene_to_family != UNASSIGNED))
    print(f"{len(gene_ids)} genes x {len(axes)} axes, {n_families} families; "
          f"{n_assigned} genes assigned, {len(gene_ids) - n_assigned} unassigned")
    if n_missing:
        print(f"note: {n_missing} family members have no expression row and were skipped")

    if args.orthologs is None:
        centroids = group_centroid_all(expression_vectors, gene_to_family, n_families)
        contributing = gene_to_family != UNASSIGNED
        print("centroids taken over ALL genes of each family")
    else:
        ortholog_mask = read_ortholog_mask(args.orthologs, gene_ids)
        centroids = group_centroid_orthologs(
            expression_vectors, gene_to_family, n_families, ortholog_mask
        )
        contributing = (gene_to_family != UNASSIGNED) & ortholog_mask
        print(f"centroids taken over ORTHOLOGS only "
              f"({int(np.count_nonzero(contributing))} contributing genes)")

    # A family that nothing contributed to gets a zero vector, not an error and
    # not a flag -- and a zero centroid is a real point at the origin, so every
    # distance later measured to it is just the gene's own norm. Report it.
    n_contributors = np.bincount(
        gene_to_family[contributing], minlength=n_families + 1
    )[1:]
    empty = np.flatnonzero(n_contributors == 0)
    if empty.size:
        print(f"WARNING: {empty.size} of {n_families} families have no contributing gene "
              f"and received a zero centroid, e.g. "
              + ", ".join(family_ids[i] for i in empty[:3]))
    else:
        print(f"every family has at least one contributing gene "
              f"(smallest: {n_contributors.min()}, largest: {n_contributors.max()})")

    write_labeled_matrix(args.out, "family_id", family_ids, axes, centroids)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
