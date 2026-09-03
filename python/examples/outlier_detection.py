#!/usr/bin/env python3
"""Find the genes whose expression has diverged from their gene family.

The chain is four calls:

    group_centroid_all          a reference per family
    distance_to_centroid        how far each gene sits from its own family's
    detect_outliers             which of those distances are extreme
    compute_shift_vector_field  the direction each divergence points in

``detect_outliers`` does not compare a gene's distance against a fixed cutoff.
It divides each distance by its family's own spread -- giving the Relative
Distance Index, so that a wide family and a tight one are on the same footing
-- and then flags the top ``--percentile`` of those indices. The scale of a
family whose own spread is too small to estimate comes from a LOESS fit of
spread against mean distance over all families.

Genes belonging to no family are carried through the whole chain: they are
never flagged and never enter the threshold.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/outlier_detection.py

which uses the bundled 458 families over ``material/normalization.tsv`` and
writes ``results/outliers.tsv`` plus the shift vectors of the flagged genes.

Reading the input files is plumbing and lives in ``example_data.py``; what is
below is the workflow itself.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from example_data import (UNASSIGNED, build_gene_to_family, read_expression,
                          read_families, read_ortholog_mask, write_labeled_matrix)
from tensor_omics import (compute_shift_vector_field, detect_outliers,
                          distance_to_centroid, group_centroid_all,
                          group_centroid_orthologs)

DEFAULT_FAMILIES = Path("material/filtered_families.tsv")
DEFAULT_EXPRESSION = Path("material/normalization.tsv")
DEFAULT_OUTPUT = Path("results/outliers.tsv")


def write_outlier_table(out_path, gene_ids, family_ids, gene_to_family,
                        distances, quantile, is_outlier):
    """Write one row per gene that belongs to a family."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write("gene_id\tfamily_id\tdistance_to_centroid\tquantile\tis_outlier\n")
        for i_gene, gene_id in enumerate(gene_ids):
            i_family = gene_to_family[i_gene]
            if i_family == UNASSIGNED:
                continue
            handle.write(
                f"{gene_id}\t{family_ids[i_family - 1]}\t{distances[i_gene]:.6g}\t"
                f"{quantile[i_gene]:.6g}\t{int(is_outlier[i_gene])}\n"
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES,
                        help=f"family table (default: {DEFAULT_FAMILIES})")
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION,
                        help=f"normalized gene x axis table (default: {DEFAULT_EXPRESSION})")
    parser.add_argument("--orthologs", type=Path, default=None,
                        help="gene ids of the orthologs, one per line; if given, the "
                             "centroids are taken over these genes only")
    parser.add_argument("--percentile", type=float, default=0.95,
                        help="fraction, NOT a percentage: 0.95 flags the top 5%% (default: 0.95)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT,
                        help=f"where to write the per-gene table (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args(argv)

    family_ids, genes_of_family = read_families(args.families)
    gene_ids, axes, expression_vectors = read_expression(args.expression)
    gene_to_family, _ = build_gene_to_family(gene_ids, genes_of_family)
    n_families = len(family_ids)

    assigned = gene_to_family != UNASSIGNED
    print(f"{len(gene_ids)} genes x {len(axes)} axes; "
          f"{int(assigned.sum())} in {n_families} families")

    if args.orthologs is None:
        centroids = group_centroid_all(expression_vectors, gene_to_family, n_families)
    else:
        centroids = group_centroid_orthologs(
            expression_vectors, gene_to_family, n_families,
            read_ortholog_mask(args.orthologs, gene_ids),
        )

    # Genes with no family get the -1 distance sentinel rather than a number.
    distances = distance_to_centroid(expression_vectors, centroids, gene_to_family)
    real = distances[assigned]
    print(f"distance to own centroid: median {np.median(real):.4g}, max {real.max():.4g}")

    result = detect_outliers(n_families, distances, gene_to_family,
                             percentile=args.percentile)
    is_outlier, quantile = result["is_outlier"], result["quantile"]

    n_flagged = int(is_outlier.sum())
    print(f"flagged {n_flagged} of {int(assigned.sum())} assigned genes "
          f"({n_flagged / max(int(assigned.sum()), 1):.2%}) at percentile {args.percentile}")
    if is_outlier[~assigned].any():
        print("WARNING: a gene with no family was flagged; this should not happen")

    # The field holds two vectors per gene: its family's centroid, then the shift
    # from that centroid to the gene. The shift is what carries the direction of
    # the divergence, which the distance alone throws away.
    field = compute_shift_vector_field(expression_vectors, centroids, gene_to_family)
    shift_vectors = field[:, 1, :]

    write_outlier_table(args.out, gene_ids, family_ids, gene_to_family,
                        distances, quantile, is_outlier)
    print(f"wrote {args.out}")

    flagged = np.flatnonzero(is_outlier)
    shift_path = args.out.with_name(args.out.stem + "_shift_vectors.tsv")
    write_labeled_matrix(shift_path, "gene_id", [gene_ids[i] for i in flagged],
                         axes, shift_vectors[:, flagged])
    print(f"wrote {shift_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
