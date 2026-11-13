#!/usr/bin/env python3
"""Turn OrthoFinder's output into the family inputs Tensor Omics needs.

OrthoFinder itself runs outside Tensor Omics. What this script covers is the
handoff: reading its results, checking them against your expression table, and
deriving the one array TOX needs but no bundled file provides -- the ortholog
mask that `group_centroid_orthologs` takes.

Two OrthoFinder outputs matter here, and they answer different questions:

  Orthogroups/Orthogroups.tsv   the FAMILIES: one row per orthogroup, one
                               column per species, genes listed per cell
  Orthologues/<sp>/<a>__v__<b>.tsv   the ORTHOLOGY relations: three columns,
                               Orthogroup and one per species, each row the
                               gene(s) in one species orthologous to the
                               gene(s) in the other

A row with exactly one gene on each side is a one-to-one orthology. Genes that
appear only in one-to-many or many-to-many rows are the duplicated copies --
the ones an SND analysis is about. Treating the one-to-one participants as
"the orthologs" is a choice, and this script implements that one; see the
recipe for what else is defensible.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/gene_family_detection.py

which reports the bundled ``material/Orthogroups.tsv`` against
``material/kallisto_sex_data_no_na.tsv``. Add ``--orthologues DIR`` to also
write the ortholog list.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from tensor_omics import (get_unassigned_mask, read_gene_ids_from_tsv_file,
                          read_orthofinder_file, validate_gene_to_family_mapping,
                          validate_string_array_uniqueness)

DEFAULT_FAMILIES = Path("material/Orthogroups.tsv")
DEFAULT_EXPRESSION = Path("material/kallisto_sex_data_no_na.tsv")


def survey_orthogroups(tsv_path):
    """Report the shape of an Orthogroups.tsv without interpreting it.

    Returns
    -------
    species : list of str
    sizes : list of int
        Genes per orthogroup, in file order.
    n_species_present : list of int
        How many species contribute at least one gene, per orthogroup.
    """
    with tsv_path.open() as handle:
        species = handle.readline().rstrip("\n").split("\t")[1:]
        sizes, present = [], []
        for line in handle:
            cells = line.rstrip("\n").split("\t")[1:]
            per_species = [[g for g in cell.split(",") if g.strip()] for cell in cells]
            sizes.append(sum(len(g) for g in per_species))
            present.append(sum(1 for g in per_species if g))
    return species, sizes, present


def one_to_one_orthologs(orthologues_dir):
    """Collect the genes taking part in at least one 1:1 orthology relation.

    Parameters
    ----------
    orthologues_dir : Path
        OrthoFinder's ``Orthologues/`` directory: one sub-directory per species,
        each holding ``<a>__v__<b>.tsv`` files with three columns -- Orthogroup,
        then one column of gene(s) per species.

    Returns
    -------
    set of str
        Genes that appear alone on both sides of at least one row. A gene in
        only one-to-many or many-to-many rows is excluded: it is one of several
        copies, which is exactly what the divergence analyses test.
    """
    orthologs = set()
    for tsv_path in sorted(orthologues_dir.rglob("*__v__*.tsv")):
        with tsv_path.open() as handle:
            next(handle)  # Orthogroup, species A, species B
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                left = [g.strip() for g in fields[1].split(",") if g.strip()]
                right = [g.strip() for g in fields[2].split(",") if g.strip()]
                if len(left) == 1 and len(right) == 1:
                    orthologs.update(left + right)
    return orthologs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--families", type=Path, default=DEFAULT_FAMILIES,
                        help=f"OrthoFinder Orthogroups.tsv (default: {DEFAULT_FAMILIES})")
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION,
                        help=f"expression table to check against (default: {DEFAULT_EXPRESSION})")
    parser.add_argument("--orthologues", type=Path, default=None,
                        help="OrthoFinder Orthologues/ directory; if given, the 1:1 "
                             "ortholog genes are written to --ortholog-out")
    parser.add_argument("--ortholog-out", type=Path, default=Path("results/orthologs.txt"),
                        help="where to write the ortholog gene list")
    args = parser.parse_args(argv)

    species, sizes, present = survey_orthogroups(args.families)
    n_families = len(sizes)
    print(f"{n_families} orthogroups over {len(species)} species: "
          + ", ".join(species))
    print(f"  family size: min {min(sizes)}, median {int(np.median(sizes))}, max {max(sizes)}")
    complete = sum(1 for p in present if p == len(species))
    print(f"  {complete} families have all {len(species)} species, "
          f"{n_families - complete} are missing at least one")
    print(f"  {sum(1 for s in sizes if s > len(species))} families hold more genes than "
          f"species, i.e. contain duplications")

    # --- the handoff: does the family table match the expression table? ------
    n_genes = sum(1 for _ in args.expression.open()) - 1
    gene_id_width = max(len(line.split("\t", 1)[0])
                        for line in args.expression.open().readlines()[1:])
    gene_ids = read_gene_ids_from_tsv_file(
        str(args.expression), gene_id_width, n_genes, n_header_rows=1, gene_col=1)

    families = read_orthofinder_file(
        str(args.families), gene_ids,
        max(len(f) for f in [line.split("\t", 1)[0] for line in args.families.open()][1:]),
        n_families, n_genes)
    gene_to_fam = np.asarray(families["gene_to_fam"])

    validate_string_array_uniqueness(gene_ids)
    validate_string_array_uniqueness(families["family_ids"])
    validate_gene_to_family_mapping(gene_to_fam, n_families)

    n_assigned = get_unassigned_mask(gene_to_fam)["n_genes_kept"]
    print(f"\nagainst {args.expression.name}: {n_assigned} of {n_genes} genes "
          f"land in a family, {n_genes - n_assigned} carry the unassigned sentinel")
    empty = n_families - len(set(gene_to_fam[gene_to_fam > 0]))
    print(f"  {empty} of {n_families} families have no gene in the expression table")

    # --- the ortholog mask, which no TOX routine produces --------------------
    if args.orthologues is None:
        print("\nno --orthologues given, so no ortholog list was written; "
              "group_centroid_orthologs needs one")
        return 0

    orthologs = one_to_one_orthologs(args.orthologues)
    in_expression = [g for g in gene_ids if g in orthologs]
    args.ortholog_out.parent.mkdir(parents=True, exist_ok=True)
    args.ortholog_out.write_text("\n".join(in_expression) + "\n")
    print(f"\n{len(orthologs)} genes take part in a 1:1 orthology; "
          f"{len(in_expression)} of them are in the expression table")
    print(f"wrote {args.ortholog_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
