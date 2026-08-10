#!/usr/bin/env python3
"""
Runs run_stc.py twice on the same input: once with the parameters you give it (the
"original" run), once with estimate_stc_parameters' own proposal for k_min/k_density/
chordal_dist_max_as_prcnt_of_range/g_max/d_max actually applied as the real run parameters
(the "estimated" run) -- not just reported, unlike --estimate-parameters on its own.
`density_quantile` has no direct CLI equivalent (see run_stc.py's own params-table "no
equivalent" note) and is skipped; reconciliation_mode/min_jsi/rmse_change_max are held fixed
across both runs, since the estimator does not touch them.

Run: python3 run_stc_pair.py <input.csv> [run_stc.py options for the "original" run]
Writes <out-prefix>_original_* and <out-prefix>_estimated_* (defaults to
results/data/<stem>_original / _estimated next to this script, same as run_stc.py's own
default), plus the matching PDF/3D reports via plot_stc.R.
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Parameters estimate_stc_parameters proposes that map directly onto a run_stc.py CLI flag.
# density_quantile deliberately excluded -- no direct equivalent, see run_stc.py's own note.
ESTIMATE_TO_CLI_FLAG = {
    "estimated_k_min": "--k-min",
    "estimated_k_density": "--k-density",
    "estimated_chordal_dist_max_as_prcnt_of_range": "--chordal-dist-max-as-prcnt-of-range",
    "estimated_g_max": "--g-max",
    "estimated_d_max": "--d-max",
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_csv")
    p.add_argument("--out-prefix", default=None,
                    help="Base prefix; '_original'/'_estimated' are appended (default: "
                         "results/data/<input file stem> next to this script)")
    # Everything else is passed straight through to the "original" run_stc.py invocation.
    args, passthrough = p.parse_known_args()

    stem = os.path.splitext(os.path.basename(args.input_csv))[0]
    base_prefix = args.out_prefix or os.path.join(HERE, "results", "data", stem)
    original_prefix = f"{base_prefix}_original"
    estimated_prefix = f"{base_prefix}_estimated"
    run_stc_py = os.path.join(HERE, "run_stc.py")

    print(f"[run_stc_pair] {stem}: running the original parameters...")
    subprocess.run([sys.executable, run_stc_py, args.input_csv, "--out-prefix", original_prefix,
                    "--estimate-parameters", *passthrough], check=True)

    with open(f"{original_prefix}_params.json") as fh:
        original_params = json.load(fh)

    estimated_flags = []
    for key, flag in ESTIMATE_TO_CLI_FLAG.items():
        if key in original_params:
            estimated_flags.extend([flag, str(original_params[key])])
    # Hold reconciliation behavior (and RMSE_change_max, which estimate_stc_parameters does
    # not propose a value for, see misc/mod_STC.md) fixed across both runs -- the estimator
    # does not touch these, so comparing the two runs makes sense only if they match.
    for flag in ("--reconciliation-mode", "--min-jsi", "--max-group-size", "--rmse-change-max"):
        key = flag.lstrip("-").replace("-", "_")
        if original_params.get(key) not in (None, "None"):
            estimated_flags.extend([flag, str(original_params[key])])

    print(f"[run_stc_pair] {stem}: running with estimate_stc_parameters' own proposal -- "
          f"{' '.join(estimated_flags)}")
    subprocess.run([sys.executable, run_stc_py, args.input_csv, "--out-prefix", estimated_prefix,
                    *estimated_flags], check=True)

    for prefix in (original_prefix, estimated_prefix):
        print(f"[run_stc_pair] {stem}: plotting {prefix}...")
        subprocess.run(["Rscript", os.path.join(HERE, "plot_stc.R"), prefix], check=True)
        subprocess.run([sys.executable, os.path.join(HERE, "render_interactive.py"), prefix], check=True)

    print(f"[run_stc_pair] {stem}: done -- {original_prefix}.pdf vs {estimated_prefix}.pdf "
          f"(and their _interactive.html siblings)")


if __name__ == "__main__":
    main()
