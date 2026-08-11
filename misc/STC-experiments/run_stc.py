#!/usr/bin/env python3
"""
Runs the STC pipeline (seeds -> ensemble_identification_merged -> ensemble_reconciliation)
on one input point-cloud CSV and writes a set of relational output CSVs that
`plot_stc.R` reads to build a PDF report.

This is the STC counterpart of branch `smoothing`'s `test_aux/test_lomanle.f90` +
`run_lomanle_tests.sh` combination -- except there is no need for a standalone compiled
Fortran driver here: STC's Python bindings (`python/tensor_omics`) already expose the whole
pipeline, so this script calls them directly instead of shelling out to a Fortran binary.

Run: python3 run_stc.py <input.csv> [options]
See `--help` for the full option list, or misc/STC-experiments/README.md for a walkthrough.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO_ROOT, "python"))

from tensor_omics import (build_kd_index, seeds, ensemble_identification_merged,  # noqa: E402
                          ensemble_reconciliation, estimate_stc_parameters)
from tensor_omics.error_handling import ToxError  # noqa: E402

STOP_REASON_NAMES = {
    0: "error",
    1: "max_size",
    2: "rejected_after_stable",
    3: "rejected_immediately",
    4: "fixed_point",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_csv", help="Point-cloud CSV; every numeric column is treated as one ambient dimension")
    p.add_argument("--out-prefix", default=None,
                    help="Prefix for the output CSVs (default: results/<input file stem> next to this script)")
    p.add_argument("--k", type=int, default=None,
                    help="Shortcut: sets both --k-min and --k-density at once (individual flags still override)")
    p.add_argument("--k-min", type=int, default=None,
                    help="Growth radius: neighborhood size for the median-distance radius (default 30, or --k)")
    p.add_argument("--k-density", type=int, default=None,
                    help="Seeding: neighborhood size for both the density estimate and the coverage radius, "
                         "see density_labels/calc_ensemble_growth_radius (default 30, or --k -- the same "
                         "default value k_min has, not a live reference: seeds() and "
                         "ensemble_identification_merged() are separate top-level calls, so --k is the only "
                         "way to move both together from the command line)")
    p.add_argument("--exclusion-radius-percentile", type=float, default=50.0,
                    help="Seeding: percentile (0-100) of the k_density neighbor distances used as each seed's "
                         "coverage/exclusion radius, see seeds (default 50.0, the median)")
    p.add_argument("--bandwidth-percentile", type=float, default=68.27,
                    help="Seeding: percentile (0-100) of the k_density neighbor distances used as the local "
                         "Gaussian KDE bandwidth for density ranking, see density_labels (default 68.27)")
    p.add_argument("--chordal-dist-max-as-prcnt-of-range", type=float, default=0.5,
                    help="Accept: maximum tolerated chordal distance between tangent bases, as a fraction "
                         "(0-1) of its own [0, sqrt(d)] range -- see accept_ensemble's tangent-space-drift "
                         "criterion in misc/mod_STC.md (default 0.5, ~ sin(30deg), a literal carry-over of "
                         "the old single-angle default under the d=1 special case)")
    p.add_argument("--d-max", type=int, default=1,
                    help="Accept: maximum tolerated change in intrinsic dimension between growth steps (default 1)")
    p.add_argument("--g-max", type=float, default=3.0,
                    help="Accept: maximum tolerated |log(G_tp1/G_t)| in the spectral gap (default 3.0)")
    p.add_argument("--rmse-change-max", type=float, default=abs(np.log(1.5)),
                    help="Accept: maximum tolerated |log(RMSE_tp1/RMSE_t)| in the residual, see "
                         "accept_ensemble's residual-drift criterion (default |log(1.5)| ~ 0.405)")
    p.add_argument("--f-max", type=float, default=0.95,
                    help="Stop condition 1: ensemble size fraction of N above which growth is abandoned (default 0.95)")
    p.add_argument("--a", type=int, default=2,
                    help="Stop condition 2: accepted-iteration count for a later rejection to count as stable (default 2)")
    p.add_argument("--o", type=int, default=10,
                    help="Trailing observable-history window depth (default 10)")
    p.add_argument("--reconciliation-mode", choices=["report", "merge_overlap_coefficient", "merge_any"],
                    default="merge_overlap_coefficient",
                    help="How Ensemble Reconciliation processes intersections (default merge_overlap_coefficient)")
    p.add_argument("--min-overlap-coefficient", type=float, default=0.9,
                    help="Reconciliation: minimum Overlap Coefficient (|intersection| / min(|A|,|B|)) for mode "
                         "merge_overlap_coefficient (default 0.9)")
    p.add_argument("--max-group-size", type=int, default=None,
                    help="Reconciliation: max ensembles per super-ensemble (default min(1024, n_ensembles))")
    p.add_argument("--estimate-parameters", action="store_true",
                    help="Also run estimate_stc_parameters and report its k_min/k_density/density_quantile/"
                         "chordal_dist_max_as_prcnt_of_range/G_max/d_max estimates in params.json/.txt -- "
                         "purely informational, does not change which parameters this run itself actually "
                         "uses (see misc/mod_STC.md, 'Estimate parameters from data'). RMSE_change_max is "
                         "not estimated by this SKG -- not part of the reported estimate.")
    p.add_argument("--n-anchors", type=int, default=5,
                    help="Parameter estimation: number of estimator anchors, see estimate_stc_parameters "
                         "(default 5, only used with --estimate-parameters)")
    p.add_argument("--seed-max-set-size", type=float, default=5.0,
                    help="Parameter estimation: percent of N at which estimator-anchor growth stops, see "
                         "estimate_stc_parameters (default 5.0, only used with --estimate-parameters)")
    p.add_argument("--first-quartile-percentile", type=float, default=25.0,
                    help="Parameter estimation: percentile of the pairwise-EA-comparison distributions used "
                         "for chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters "
                         "(default 25.0, only used with --estimate-parameters)")
    args = p.parse_args()
    # --k-min/--k-density, in that order, always win over --k; --k only fills in whichever of
    # the two was not given explicitly; 30 is the hardcoded fallback if neither was.
    args.k_min = args.k_min if args.k_min is not None else (args.k if args.k is not None else 30)
    args.k_density = args.k_density if args.k_density is not None else (args.k if args.k is not None else 30)
    return args


def read_points(csv_path):
    """Every numeric column of the CSV is one ambient dimension -- matches how
    r/plot_lomanle_spheres.R on branch `smoothing` folds x1, x2, y_original into x, y[, z]:
    STC does not distinguish a "response" column from a "feature" column, it just clusters
    points in ambient space, so there is no reason to treat them differently here either."""
    df = pd.read_csv(csv_path)
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if len(numeric_cols) < 2:
        raise ValueError(f"{csv_path}: need at least 2 numeric columns, found {numeric_cols}")
    points = df[numeric_cols].to_numpy(dtype=np.float64)
    return points, numeric_cols


def last_valid_history_index(k_history_col):
    """Index of the rightmost nonzero entry in one ensemble's k_history column, or None if
    the ensemble never had an observable computed at all (e.g. a size-1 trivial ensemble --
    see ensemble_identification_kernel's own "isolated seed" note)."""
    nonzero = np.nonzero(k_history_col)[0]
    return int(nonzero[-1]) if len(nonzero) > 0 else None


def main():
    args = parse_args()

    stem = os.path.splitext(os.path.basename(args.input_csv))[0]
    out_prefix = args.out_prefix or os.path.join(HERE, "results", "data", stem)
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    points, dim_names = read_points(args.input_csv)
    n_vectors, n_dimensions = points.shape
    vectors = np.asfortranarray(points.T)  # (n_dimensions, n_vectors), as the bindings expect

    dimension_order = np.arange(1, n_dimensions + 1, dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)

    # --- optional: estimate_stc_parameters, purely informational (see --estimate-parameters
    # help text above) -- reported in params.json/.txt, never substituted for this run's own
    # k_min/k_density/etc.
    estimated_params = {}
    if args.estimate_parameters:
        try:
            est = estimate_stc_parameters(
                vectors, kd_indices, dimension_order,
                n_anchors=args.n_anchors, seed_max_set_size=args.seed_max_set_size,
                first_quartile_percentile=args.first_quartile_percentile,
            )
            estimated_params = {
                "estimated_k_min": round(float(est["estimated_k_min"])),
                "estimated_k_density": round(float(est["estimated_k_density"])),
                "estimated_density_quantile": float(est["estimated_density_quantile"]),
                "estimated_chordal_dist_max_as_prcnt_of_range":
                    float(est["estimated_chordal_dist_max_as_prcnt_of_range"]),
                "estimated_g_max": float(est["estimated_G_max"]),
                "estimated_d_max": round(float(est["estimated_d_max"])),
            }
            print(f"[run_stc] {stem}: estimate_stc_parameters -> "
                  f"k_min={estimated_params['estimated_k_min']}, "
                  f"k_density={estimated_params['estimated_k_density']}, "
                  f"density_quantile={estimated_params['estimated_density_quantile']:.4g}, "
                  f"chordal_dist_max_as_prcnt_of_range="
                  f"{estimated_params['estimated_chordal_dist_max_as_prcnt_of_range']:.3g}, "
                  f"g_max={estimated_params['estimated_g_max']:.4g}, "
                  f"d_max={estimated_params['estimated_d_max']}")
        except ToxError as error:
            print(f"[run_stc] {stem}: estimate_stc_parameters failed ({error}); "
                  f"continuing without an estimate (see misc/mod_STC.md's own note that this is a "
                  f"heuristic that can genuinely fail on some inputs, e.g. too few usable anchors)")

    seed_selection_mask = seeds(vectors, kd_indices, dimension_order, k_density=args.k_density,
                                exclusion_radius_percentile=args.exclusion_radius_percentile,
                                bandwidth_percentile=args.bandwidth_percentile)
    n_ensembles = int(np.count_nonzero(seed_selection_mask))
    print(f"[run_stc] {stem}: N={n_vectors}, D={n_dimensions}, seeds found={n_ensembles}")

    result = ensemble_identification_merged(
        vectors, kd_indices, dimension_order, seed_selection_mask,
        args.chordal_dist_max_as_prcnt_of_range, args.d_max, args.g_max, args.rmse_change_max, args.o,
        k_min=args.k_min, f_max=args.f_max, a=args.a,
    )

    # --- points.csv: one row per input point, wide ---
    points_df = pd.DataFrame(points, columns=dim_names)
    points_df.insert(0, "point_id", np.arange(1, n_vectors + 1))
    points_df["n_ensembles"] = result["ensemble_masks"].sum(axis=1)
    points_df["n_low_confidence_ensembles"] = result["ensemble_low_confidence_masks"].sum(axis=1)
    points_df.to_csv(f"{out_prefix}_points.csv", index=False)

    # --- low_confidence_membership.csv: long, one row per (point, ensemble) whose iteration-1
    # bootstrap mask covers it -- reported for every seed regardless of stop_reason, see
    # low_confidence_mask/ensemble_low_confidence_masks in misc/mod_STC.md, "Ensemble
    # identification", "Output". A point with n_ensembles==0 but that appears here has a
    # low-confidence fallback available; one that appears in neither has none at all.
    low_confidence_rows = []
    for e in range(n_ensembles):
        member_points = np.nonzero(result["ensemble_low_confidence_masks"][:, e])[0] + 1
        for pid in member_points:
            low_confidence_rows.append((pid, e + 1))
    low_confidence_df = pd.DataFrame(low_confidence_rows, columns=["point_id", "ensemble_id"])
    low_confidence_df.to_csv(f"{out_prefix}_low_confidence_membership.csv", index=False)

    # --- membership.csv: long, one row per (point, ensemble) that actually contains it ---
    seed_indices = np.nonzero(seed_selection_mask)[0] + 1  # 1-indexed, matches ensemble_id below
    membership_rows = []
    for e in range(n_ensembles):
        member_points = np.nonzero(result["ensemble_masks"][:, e])[0] + 1
        for pid in member_points:
            membership_rows.append((pid, e + 1, int(pid == seed_indices[e])))
    membership_df = pd.DataFrame(membership_rows, columns=["point_id", "ensemble_id", "is_seed"])
    membership_df.to_csv(f"{out_prefix}_membership.csv", index=False)

    # --- ensembles.csv: one row per ensemble, its final retained state ---
    ensemble_rows = []
    for e in range(n_ensembles):
        last = last_valid_history_index(result["ensemble_k_history"][:, e])
        row = {
            "ensemble_id": e + 1,
            "seed_point_id": int(seed_indices[e]),
            "stop_reason": STOP_REASON_NAMES.get(int(result["ensemble_stop_reason"][e]), "unknown"),
            "growth_radius": float(result["ensemble_growth_radii"][e]),
            "size": int(result["ensemble_masks"][:, e].sum()),
        }
        for j, name in enumerate(dim_names):
            row[f"mu_{name}"] = np.nan
        row["d"] = 0
        row["G"] = np.nan
        for t in (1, 2):
            for name in dim_names:
                row[f"u{t}_{name}"] = 0.0
            row[f"s{t}"] = 0.0

        if last is not None:
            d = int(result["ensemble_d_history"][last, e])
            k = int(result["ensemble_k_history"][last, e])
            mu = result["ensemble_mu_history"][:, last, e]
            S = result["ensemble_S_history"][:, last, e]
            U = result["ensemble_U_history"][:, :, last, e]
            row["d"] = d
            row["G"] = float(result["ensemble_G_history"][last, e])
            for j, name in enumerate(dim_names):
                row[f"mu_{name}"] = float(mu[j])
            # tangent_scale_j = sqrt(eigenvalue_j) = S_j / sqrt(k - 1), see
            # tox_shape_truthful_clustering_observable_kernel's own normal_error/tangent_scales
            for t in range(1, min(d, 2) + 1):
                scale = float(S[t - 1] / np.sqrt(max(k - 1, 1)))
                row[f"s{t}"] = scale
                for j, name in enumerate(dim_names):
                    row[f"u{t}_{name}"] = float(U[j, t - 1])
        ensemble_rows.append(row)
    ensembles_df = pd.DataFrame(ensemble_rows)
    ensembles_df.to_csv(f"{out_prefix}_ensembles.csv", index=False)

    # --- reconciliation: super_ensembles.csv / super_ensembles_overlap_coefficient.csv ---
    group_rows, overlap_coefficient_rows = [], []
    if n_ensembles >= 2:
        max_group_size = args.max_group_size or min(1024, n_ensembles)
        recon = ensemble_reconciliation(
            result["ensemble_masks"], max_group_size,
            mode=args.reconciliation_mode, min_overlap_coefficient=args.min_overlap_coefficient,
            report_overlap_coefficient=True,
        )
        n_groups = int(recon["n_super_ensembles"])
        super_ensembles = recon["super_ensembles"]
        super_ensembles_overlap_coefficient = recon["super_ensembles_overlap_coefficient"]
        for g in range(n_groups):
            members = super_ensembles[:, g]
            members = members[members > 0]
            for m in members:
                group_rows.append((g + 1, int(m)))
            for row_idx in range(len(members) - 1):
                overlap_coefficient_rows.append((g + 1, int(members[row_idx]), int(members[row_idx + 1]),
                                  float(super_ensembles_overlap_coefficient[row_idx, g])))
        print(f"[run_stc] {stem}: reconciliation mode={args.reconciliation_mode} -> {n_groups} super-ensemble(s)")
    else:
        print(f"[run_stc] {stem}: fewer than 2 ensembles, skipping reconciliation")

    pd.DataFrame(group_rows, columns=["group_id", "ensemble_id"]).to_csv(f"{out_prefix}_super_ensembles.csv", index=False)
    pd.DataFrame(overlap_coefficient_rows, columns=["group_id", "ensemble_id_from", "ensemble_id_to", "overlap_coefficient"]).to_csv(
        f"{out_prefix}_super_ensembles_overlap_coefficient.csv", index=False)

    # --- ensemble_overlap_coefficient_matrix.csv: every pairwise Overlap Coefficient between
    # non-empty ensembles, not just the consecutive-within-a-group pairs
    # super_ensembles_overlap_coefficient.csv reports -- for the "Overlap Coefficient between
    # ensembles" heatmap, which needs the full matrix (including pairs that never qualified
    # for reconciliation at all, OC=0) rather than only super-ensemble members. Cheap: a
    # single boolean matmul over the same ensemble_masks reconciliation already used.
    nonempty = np.nonzero(result["ensemble_masks"].sum(axis=0) > 0)[0]
    overlap_coefficient_matrix_rows = []
    if len(nonempty) >= 2:
        masks_i = result["ensemble_masks"][:, nonempty].astype(np.int64)
        sizes_i = masks_i.sum(axis=0)
        intersection = masks_i.T @ masks_i
        min_size = np.minimum(sizes_i[:, None], sizes_i[None, :])
        oc_full = np.divide(intersection, min_size, out=np.zeros_like(intersection, dtype=np.float64), where=min_size > 0)
        for a_idx, e_a in enumerate(nonempty):
            for b_idx, e_b in enumerate(nonempty):
                if e_b <= e_a:
                    continue
                overlap_coefficient_matrix_rows.append((int(e_a) + 1, int(e_b) + 1, float(oc_full[a_idx, b_idx])))
    pd.DataFrame(overlap_coefficient_matrix_rows, columns=["ensemble_id_1", "ensemble_id_2", "overlap_coefficient"]).to_csv(
        f"{out_prefix}_ensemble_overlap_coefficient_matrix.csv", index=False)

    # --- params.json / params.txt: for the plot script's caption. Both are written: .json
    # for anything that wants to parse it programmatically, .txt (plain key=value lines) so
    # plot_stc.R does not need the jsonlite package, which is not installed in every R
    # environment this might run in (this one included).
    params = vars(args).copy()
    params["n_vectors"] = n_vectors
    params["n_dimensions"] = n_dimensions
    params["n_ensembles"] = n_ensembles
    params["dim_names"] = dim_names
    params.update(estimated_params)
    with open(f"{out_prefix}_params.json", "w") as fh:
        json.dump(params, fh, indent=2)
    with open(f"{out_prefix}_params.txt", "w") as fh:
        for key, value in params.items():
            if isinstance(value, list):
                value = ",".join(str(v) for v in value)
            fh.write(f"{key}={value}\n")

    print(f"[run_stc] wrote {out_prefix}_{{points,membership,low_confidence_membership,ensembles,"
          f"super_ensembles,super_ensembles_overlap_coefficient,ensemble_overlap_coefficient_matrix,"
          f"params}}.{{csv,json,txt}}")


if __name__ == "__main__":
    main()
