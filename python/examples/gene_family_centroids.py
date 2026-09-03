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

This script is deliberately self-contained: copy it and adapt it.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from tensor_omics import group_centroid_all, group_centroid_orthologs

#: gene_to_family entry for a gene that belongs to no family (M_GENE_TO_FAM_SENTINEL)
UNASSIGNED = 0

DEFAULT_FAMILIES = Path("material/filtered_families.tsv")
DEFAULT_EXPRESSION = Path("material/normalization.tsv")
DEFAULT_OUTPUT = Path("results/family_centroids.tsv")


def read_families(tsv_path):
    """Read an OrthoFinder-style family table.

    Parameters
    ----------
    tsv_path : Path
        First column is the family id; every further column holds that family's
        gene ids for one species, comma-separated and possibly empty.

    Returns
    -------
    family_ids : list of str
        Family names, in file order. Their positions define the family indices
        used below, which are 1-based because Fortran indexes from 1.
    genes_of_family : list of list of str
        Member gene ids per family, aligned with `family_ids`.
    """
    family_ids, genes_of_family = [], []
    with tsv_path.open() as handle:
        next(handle)  # header: family id + one column per species
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            members = [
                gene.strip()
                for cell in fields[1:]
                for gene in cell.split(",")
                if gene.strip()
            ]
            family_ids.append(fields[0])
            genes_of_family.append(members)
    return family_ids, genes_of_family


def read_expression(tsv_path):
    """Read a gene x axis expression table into TOX's (n_axes, n_genes) layout.

    Returns
    -------
    gene_ids : list of str
    axes : list of str
        Column names, i.e. the tissues/conditions spanning the space.
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes), order='F'
        One gene per column, because Fortran is column-major and TOX stores
        vectors as columns.
    """
    with tsv_path.open() as handle:
        axes = handle.readline().rstrip("\n").split("\t")[1:]
        gene_ids, rows = [], []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene_ids.append(fields[0])
            rows.append([float(value) for value in fields[1:]])

    # rows is gene-major; TOX wants one gene per column, so transpose.
    return gene_ids, axes, np.asfortranarray(np.array(rows, dtype=np.float64).T)


def build_gene_to_family(gene_ids, family_ids, genes_of_family):
    """Map each expression row onto a 1-based family index, 0 where unassigned.

    A gene listed in several families keeps its first assignment, and a family
    member absent from the expression table is simply skipped -- both are
    reported by the caller rather than silently dropped.

    Returns
    -------
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
    n_missing : int
        Family members that have no row in the expression table.
    """
    row_of_gene = {gene: row for row, gene in enumerate(gene_ids)}
    gene_to_family = np.full(len(gene_ids), UNASSIGNED, dtype=np.int32)

    n_missing = 0
    for i_family, members in enumerate(genes_of_family, start=1):
        for gene in members:
            row = row_of_gene.get(gene)
            if row is None:
                n_missing += 1
            elif gene_to_family[row] == UNASSIGNED:
                gene_to_family[row] = i_family

    return gene_to_family, n_missing


def read_ortholog_mask(txt_path, gene_ids):
    """Read a list of ortholog gene ids into a mask over `gene_ids`."""
    orthologs = {line.strip() for line in txt_path.open() if line.strip()}
    return np.array([gene in orthologs for gene in gene_ids], dtype=bool)


def write_centroids(out_path, family_ids, axes, centroids):
    """Write the (n_axes, n_families) centroid matrix as a family-per-row TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write("family_id\t" + "\t".join(axes) + "\n")
        for i_family, family_id in enumerate(family_ids):
            values = "\t".join(f"{value:.6g}" for value in centroids[:, i_family])
            handle.write(f"{family_id}\t{values}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES,
                        help=f"family table (default: {DEFAULT_FAMILIES})")
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION,
                        help=f"normalized gene x axis table (default: {DEFAULT_EXPRESSION})")
    parser.add_argument("--orthologs", type=Path, default=None,
                        help="gene ids of the orthologs, one per line; "
                             "given, centroids are taken over these genes only")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT,
                        help=f"where to write the centroids (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args(argv)

    family_ids, genes_of_family = read_families(args.families)
    gene_ids, axes, expression_vectors = read_expression(args.expression)
    gene_to_family, n_missing = build_gene_to_family(gene_ids, family_ids, genes_of_family)

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

    write_centroids(args.out, family_ids, axes, centroids)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
