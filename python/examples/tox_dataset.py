#!/usr/bin/env python3
"""Build a TOX dataset from tabular input and store it as a .txdata archive.

Tensor Omics keeps its data as plain Fortran arrays tied together by position:
column *i* of ``expression`` is the gene named by ``gene_ids[i]``, and entry *i*
of ``gene_to_family`` names that gene's family in ``family_ids``. There is no
object holding them -- the arrays are yours to keep. A ``.txdata`` archive is
how you keep them together on disk.

This script reads the bundled tables with the library's own readers, checks the
relationships between the arrays, writes an archive, and reads it back to show
the round trip is exact.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/tox_dataset.py

which reads ``material/kallisto_sex_data_no_na.tsv`` and
``material/Orthogroups.tsv`` and writes ``results/example.txdata``.

Note the shape of step 1 below: the readers fill arrays you allocate, so every
extent -- how many genes, how wide an identifier, how many families -- has to
be known before the first read. Sizing is the caller's job.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from tensor_omics import (get_tox_data_dims, get_unassigned_mask,
                          read_expression_vectors_tsv, read_gene_ids_from_tsv_file,
                          read_orthofinder_file, read_tox_data_into, save_tox_data,
                          validate_expression_data, validate_gene_to_family_mapping,
                          validate_string_array_uniqueness)

DEFAULT_EXPRESSION = Path("material/kallisto_sex_data_no_na.tsv")
DEFAULT_FAMILIES = Path("material/Orthogroups.tsv")
DEFAULT_OUTPUT = Path("results/example.txdata")


def measure(tsv_path):
    """Count the rows and columns the readers need to be told in advance.

    Returns
    -------
    n_rows : int
        Data rows, i.e. genes.
    n_value_cols : int
        Columns after the identifier column, i.e. samples.
    id_width : int
        Longest identifier, in characters. The readers store identifiers in a
        fixed-width character array, so anything longer would be truncated.
    """
    with tsv_path.open() as handle:
        n_value_cols = len(handle.readline().rstrip("\n").split("\t")) - 1
        n_rows, id_width = 0, 0
        for line in handle:
            n_rows += 1
            id_width = max(id_width, len(line.split("\t", 1)[0]))
    return n_rows, n_value_cols, id_width


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION,
                        help=f"gene x sample TSV (default: {DEFAULT_EXPRESSION})")
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES,
                        help=f"OrthoFinder Orthogroups.tsv (default: {DEFAULT_FAMILIES})")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT,
                        help=f"archive to write (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args(argv)

    # --- 1. sizes first, because the readers fill arrays you allocate --------
    n_genes, n_samples, gene_id_width = measure(args.expression)
    n_families, _, family_id_width = measure(args.families)
    print(f"{n_genes} genes x {n_samples} samples, {n_families} families")

    # --- 2. read ------------------------------------------------------------
    gene_ids = read_gene_ids_from_tsv_file(
        str(args.expression), gene_id_width, n_genes, n_header_rows=1, gene_col=1
    )

    # Filled in place, so it must already have its final shape, dtype and
    # column-major order: one gene per column, as every TOX routine expects.
    expression = np.zeros((n_samples, n_genes), dtype=np.float64, order="F")
    read_expression_vectors_tsv(
        file_list=[str(args.expression)],
        gene_ids=gene_ids,
        expression_vectors=expression,
        n_header_rows=1,
        gene_col=1,
        value_cols=np.arange(2, n_samples + 2, dtype=np.int32),
        start_row=1,
    )

    families = read_orthofinder_file(
        str(args.families), gene_ids, family_id_width, n_families, n_genes
    )
    family_ids, gene_to_family = families["family_ids"], families["gene_to_fam"]

    n_assigned = get_unassigned_mask(gene_to_family)["n_genes_kept"]
    print(f"{n_assigned} genes assigned to a family, "
          f"{n_genes - n_assigned} carrying the unassigned sentinel")

    # --- 3. check the relationships before storing them ---------------------
    # Identifiers are primary keys: a duplicate makes the position-based links
    # between the arrays ambiguous, and nothing downstream would notice.
    validate_string_array_uniqueness(gene_ids)
    validate_string_array_uniqueness(family_ids)
    validate_gene_to_family_mapping(gene_to_family, n_families)
    validate_expression_data(expression, check_non_negative=True)
    print("gene ids unique, family ids unique, mapping in range, expression non-negative")

    # --- 4. store -----------------------------------------------------------
    # Every array needs BOTH the array and its *_file name: an array passed
    # without a name is dropped from the archive without an error.
    args.out.parent.mkdir(parents=True, exist_ok=True)
    save_tox_data(
        str(args.out),
        gene_ids=gene_ids, gene_ids_file="gene_ids.bin",
        expression=expression, expression_file="expression.bin",
        gene_to_family=gene_to_family, gene_to_family_file="gene_to_family.bin",
        family_ids=family_ids, family_ids_file="family_ids.bin",
    )
    print(f"wrote {args.out}")

    # --- 5. read it back ----------------------------------------------------
    # A member that was never written reports zero extents and comes back empty,
    # which is how an archive says "this dataset has no centroids yet".
    dims = get_tox_data_dims(str(args.out))
    print("archive holds: "
          + ", ".join(f"{k}={v}" for k, v in dims.items() if v))

    back = read_tox_data_into(str(args.out))
    exact = (
        list(back["gene_ids"]) == list(gene_ids)
        and list(back["family_ids"]) == list(family_ids)
        and np.array_equal(back["expression"], expression)
        and np.array_equal(back["gene_to_family"], gene_to_family)
    )
    print(f"round trip exact: {exact}")
    print(f"centroids not stored yet, so they read back as "
          f"{back['family_centroids'].shape}")
    return 0 if exact else 1


if __name__ == "__main__":
    sys.exit(main())
