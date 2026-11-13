#!/usr/bin/env python3
"""Normalize a replicate-level expression table into TOX expression vectors.

Runs the Tensor Omics normalization pipeline on a gene x replicate TSV and
writes the resulting gene x tissue matrix. The pipeline performs, in order:

    1. gene-wise scaling by a LOESS-stabilized standard deviation,
    2. quantile normalization across replicates   (only with --quantile),
    3. averaging of the replicates of each tissue,
    4. a log2(x + 1) transformation.

The tissue layout is derived from the header: every column named
``<tissue>_rep<n>`` belongs to ``<tissue>``, and columns of the same tissue
must be adjacent. Replicate counts may differ between tissues.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/normalization.py

which normalizes the bundled ``material/kallisto_sex_data_no_na.tsv``
(88,327 genes, 67 replicates, 13 tissues) and writes
``results/normalized_expression.tsv``.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

from tensor_omics import normalization_pipeline, read_expression_vectors_tsv, read_gene_ids_from_tsv_file

#: Trailing replicate marker in a column name, e.g. "Liver_rep3" -> "Liver"
REPLICATE_SUFFIX = re.compile(r"_rep\d+$")

DEFAULT_INPUT = Path("material/kallisto_sex_data_no_na.tsv")
DEFAULT_OUTPUT = Path("results/normalized_expression.tsv")


def read_tissue_layout(tsv_path):
    """Return the tissue names and their replicate counts, in column order.

    Parameters
    ----------
    tsv_path : Path
        Tab-separated file whose first row is a header of one gene-id column
        followed by ``<tissue>_rep<n>`` columns.

    Returns
    -------
    tissues : list of str
        Tissue names, in the order their columns appear.
    reps_per_tissue : np.ndarray[np.int32]
        Number of replicate columns per tissue, aligned with `tissues`.

    Raises
    ------
    ValueError
        If a tissue's replicate columns are not adjacent, which would make the
        replicate grouping below silently average the wrong columns together.
    """
    with tsv_path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")

    tissue_of_column = [REPLICATE_SUFFIX.sub("", name) for name in header[1:]]

    tissues, reps_per_tissue = [], []
    for tissue in tissue_of_column:
        if tissues and tissues[-1] == tissue:
            reps_per_tissue[-1] += 1
        else:
            if tissue in tissues:
                raise ValueError(
                    f"columns of tissue '{tissue}' are not adjacent in {tsv_path}; "
                    "reps_per_tissue can only describe contiguous blocks"
                )
            tissues.append(tissue)
            reps_per_tissue.append(1)

    return tissues, np.asarray(reps_per_tissue, dtype=np.int32)


def read_expression(tsv_path, n_replicates):
    """Read the whole table into a (n_replicates, n_genes) column-major array.

    Parameters
    ----------
    tsv_path : Path
        The expression table.
    n_replicates : int
        Number of value columns, i.e. all columns but the gene-id column.

    Returns
    -------
    gene_ids : list of str
        Gene identifiers, in file order.
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes), order='F'
        Expression values. TOX stores one vector per *column*, because Fortran
        is column-major; this is the layout every TOX routine expects.
    """
    n_genes = sum(1 for _ in tsv_path.open()) - 1
    gene_id_width = max(len(line.split("\t", 1)[0]) for line in tsv_path.open().readlines()[1:])

    gene_ids = read_gene_ids_from_tsv_file(
        str(tsv_path), gene_id_width, n_genes, n_header_rows=1, gene_col=1
    )

    # read_expression_vectors_tsv fills its argument in place, so it must be
    # allocated with the final shape, dtype and column-major order up front.
    expr = np.zeros((n_replicates, n_genes), dtype=np.float64, order="F")
    read_expression_vectors_tsv(
        file_list=[str(tsv_path)],
        gene_ids=gene_ids,
        expression_vectors=expr,
        n_header_rows=1,
        gene_col=1,
        value_cols=np.arange(2, n_replicates + 2, dtype=np.int32),
        start_row=1,
    )
    return gene_ids, expr


def write_matrix(out_path, gene_ids, tissues, matrix):
    """Write a (n_tissues, n_genes) matrix as a gene-per-row TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write("gene_id\t" + "\t".join(tissues) + "\n")
        for i_gene, gene_id in enumerate(gene_ids):
            values = "\t".join(f"{value:.6g}" for value in matrix[:, i_gene])
            handle.write(f"{gene_id}\t{values}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help=f"replicate-level expression TSV (default: {DEFAULT_INPUT})")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT,
                        help=f"where to write the normalized matrix (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--span", type=float, default=0.7,
                        help="LOESS span for the standard-deviation fit (default: 0.7)")
    parser.add_argument("--degree", type=int, default=2,
                        help="LOESS polynomial degree (default: 2)")
    parser.add_argument("--quantile", action="store_true",
                        help="also quantile-normalize across replicates (default: off)")
    args = parser.parse_args(argv)

    tissues, reps_per_tissue = read_tissue_layout(args.input)
    gene_ids, expr = read_expression(args.input, int(reps_per_tissue.sum()))
    print(f"read {expr.shape[1]} genes x {expr.shape[0]} replicates "
          f"in {len(tissues)} tissues: "
          + ", ".join(f"{t}({n})" for t, n in zip(tissues, reps_per_tissue)))

    normalized = normalization_pipeline(
        expr, reps_per_tissue,
        span=args.span, degree=args.degree, use_quantile=args.quantile,
    )

    # A gene whose replicates are all identical carries no mean-variance
    # information, is therefore left out of the LOESS fit, and reaches the
    # output unscaled. Report how many, so a surprising count is visible.
    n_flat = int(np.count_nonzero(expr.std(axis=0) == 0.0))
    print(f"normalized to {normalized.shape[0]} tissues x {normalized.shape[1]} genes; "
          f"value range [{normalized.min():.4g}, {normalized.max():.4g}]")
    if n_flat:
        print(f"note: {n_flat} genes have zero variance across replicates "
              "and were passed through unscaled")

    write_matrix(args.out, gene_ids, tissues, normalized)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
