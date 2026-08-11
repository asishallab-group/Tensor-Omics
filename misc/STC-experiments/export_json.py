#!/usr/bin/env python3
"""
Consolidates the CSVs run_stc.py already wrote for one dataset/parameter combination into a
single JSON file, for interactive_template.html to load. Reads only what run_stc.py already
produced -- no recomputation, no dependency on the Python bindings -- so this can run against
any prefix's output, including ones from an older run_stc.py invocation, as long as the CSVs
are there.

Run: python3 export_json.py <prefix>
  writes <prefix>.json

Points get their ensemble/low-confidence memberships, and whether they are a seed,
precomputed as lists directly on the point record -- the whole point of this export is a
single hover lookup in the browser, not another join.
"""

import csv
import json
import os
import sys
from collections import defaultdict


def read_csv_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_params_txt(path):
    params = {}
    if not os.path.exists(path):
        return params
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "=" not in line:
                continue
            key, _, value = line.partition("=")
            params[key] = value
    return params


def to_number(value):
    """Best-effort str -> int/float, leaving genuinely non-numeric strings alone -- the JSON
    export should not force e.g. reconciliation_mode or stop_reason through a numeric cast.
    Empty (pandas-NaN-written) fields become null, not a lingering empty string -- e.g. G/mu
    for an empty (size=0) ensemble, which never had an observable computed at all."""
    if value == "":
        return None
    try:
        as_float = float(value)
        as_int = int(as_float)
        return as_int if as_int == as_float else as_float
    except (TypeError, ValueError):
        return value


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 export_json.py <prefix>", file=sys.stderr)
        sys.exit(1)
    prefix = sys.argv[1]

    params = read_params_txt(f"{prefix}_params.txt")
    dim_names = params.get("dim_names", "").split(",") if params.get("dim_names") else []

    points_rows = read_csv_rows(f"{prefix}_points.csv")
    membership_rows = read_csv_rows(f"{prefix}_membership.csv")
    low_confidence_rows = read_csv_rows(f"{prefix}_low_confidence_membership.csv")
    ensembles_rows = read_csv_rows(f"{prefix}_ensembles.csv")
    super_ensembles_rows = read_csv_rows(f"{prefix}_super_ensembles.csv")
    overlap_coefficient_matrix_rows = read_csv_rows(f"{prefix}_ensemble_overlap_coefficient_matrix.csv")

    ensembles_by_point = defaultdict(list)
    seed_ensembles_by_point = defaultdict(list)
    for row in membership_rows:
        pid = int(row["point_id"])
        eid = int(row["ensemble_id"])
        ensembles_by_point[pid].append(eid)
        if row.get("is_seed") in ("1", "1.0", "True", "true"):
            seed_ensembles_by_point[pid].append(eid)

    low_confidence_by_point = defaultdict(list)
    for row in low_confidence_rows:
        low_confidence_by_point[int(row["point_id"])].append(int(row["ensemble_id"]))

    super_ensemble_by_ensemble = {}
    groups = defaultdict(list)
    for row in super_ensembles_rows:
        gid = int(row["group_id"])
        eid = int(row["ensemble_id"])
        groups[gid].append(eid)
        super_ensemble_by_ensemble[eid] = gid

    points = []
    for row in points_rows:
        pid = int(row["point_id"])
        points.append({
            "id": pid,
            "coords": [float(row[name]) for name in dim_names],
            "n_ensembles": int(float(row.get("n_ensembles", 0))),
            "n_low_confidence_ensembles": int(float(row.get("n_low_confidence_ensembles", 0))),
            "ensembles": sorted(ensembles_by_point.get(pid, [])),
            "low_confidence_ensembles": sorted(low_confidence_by_point.get(pid, [])),
            "seed_of": sorted(seed_ensembles_by_point.get(pid, [])),
        })

    ensembles = []
    for row in ensembles_rows:
        eid = int(row["ensemble_id"])
        rec = {
            "id": eid,
            "seed_point_id": int(row["seed_point_id"]),
            "stop_reason": row["stop_reason"],
            "growth_radius": float(row["growth_radius"]),
            "size": int(row["size"]),
            "d": int(row["d"]),
            "G": to_number(row["G"]),
            "mu": [to_number(row[f"mu_{name}"]) for name in dim_names],
            "super_ensemble_id": super_ensemble_by_ensemble.get(eid),
        }
        for t in (1, 2):
            s_key = f"s{t}"
            if s_key in row and float(row[s_key]) > 0:
                rec[f"u{t}"] = [to_number(row[f"u{t}_{name}"]) for name in dim_names]
                rec[s_key] = float(row[s_key])
        ensembles.append(rec)

    super_ensembles = [{"group_id": gid, "ensemble_ids": sorted(eids)} for gid, eids in sorted(groups.items())]

    overlap_coefficient_matrix = [{"a": int(row["ensemble_id_1"]), "b": int(row["ensemble_id_2"]),
                                    "overlap_coefficient": float(row["overlap_coefficient"])}
                                   for row in overlap_coefficient_matrix_rows]

    payload = {
        "dim_names": dim_names,
        "params": {k: to_number(v) for k, v in params.items()},
        "points": points,
        "ensembles": ensembles,
        "super_ensembles": super_ensembles,
        "overlap_coefficient_matrix": overlap_coefficient_matrix,
    }

    out_path = f"{prefix}.json"
    with open(out_path, "w") as fh:
        json.dump(payload, fh)
    print(f"[export_json] wrote {out_path} ({len(points)} points, {len(ensembles)} ensembles)")


if __name__ == "__main__":
    main()
