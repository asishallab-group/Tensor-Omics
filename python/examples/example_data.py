"""Reading the bundled example inputs.

This module is only plumbing: it turns the TSV files under ``material/`` into the
arrays the Tensor Omics API expects. None of it is Tensor Omics usage -- the
recipes' actual workflows live in the example scripts that import it, so that
each of those shows its analysis and nothing else.

Your own data will not be in these formats. What is worth copying from here is
the shape of the result, not the parsing:

* expression vectors are ``(n_axes, n_genes)`` and column-major, because
  Fortran is column-major and TOX stores one vector per *column*;
* ``gene_to_family`` is one 1-based family index per gene, with
  :data:`UNASSIGNED` for a gene that belongs to no family.
"""

import numpy as np

#: gene_to_family entry for a gene that belongs to no family (M_GENE_TO_FAM_SENTINEL)
UNASSIGNED = 0


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
        Family names, in file order. Their positions define the family indices,
        which are 1-based because Fortran indexes from 1.
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
        One gene per column.
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


def build_gene_to_family(gene_ids, genes_of_family):
    """Map each expression row onto a 1-based family index, UNASSIGNED where none.

    A gene listed in several families keeps its first assignment, and a family
    member absent from the expression table is skipped -- both are reported to
    the caller rather than silently dropped.

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
    """Read a list of ortholog gene ids into a boolean mask over `gene_ids`."""
    orthologs = {line.strip() for line in txt_path.open() if line.strip()}
    return np.array([gene in orthologs for gene in gene_ids], dtype=bool)


def write_labeled_matrix(out_path, label_name, labels, columns, matrix):
    """Write a (len(columns), len(labels)) matrix as a label-per-row TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        handle.write(f"{label_name}\t" + "\t".join(columns) + "\n")
        for i_label, label in enumerate(labels):
            values = "\t".join(f"{value:.6g}" for value in matrix[:, i_label])
            handle.write(f"{label}\t{values}\n")
