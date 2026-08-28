#!/bin/bash
# Smoke test for stc_sa_cli: builds the project, runs the simulated-annealing CLI end to end
# (a basic multi-dataset run, --store-iteration-shatter-results, and validation-error paths),
# and checks the expected output files/exit codes. Same convention as test_stc_cli.sh -- not
# part of `test_runner.sh`'s Fortran/Python/R suites, this exercises the C executable itself.
#
# Usage: C-layer/test_stc_sa_cli.sh [--skip-build]
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"

if [[ "${1:-}" != "--skip-build" ]]; then
  TOX_SKIP_CODE_GENERATION=1 ./build.sh
fi

BIN=$(find build -name "stc_sa_cli" -type f -executable | head -1)
if [[ -z "$BIN" ]]; then
  echo "test_stc_sa_cli.sh: could not find the built stc_sa_cli executable" >&2
  exit 1
fi

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# Two small, self-contained benchmark datasets (a noisy line, mod_STC.md's JSON schema) -- a
# suite of 2, so the multi-dataset aggregation path is exercised, not just m=1. Kept
# deliberately tiny (60 points, a 300-point reference cloud) so the whole smoke test runs in a
# few seconds; a dense-enough reference cloud only needs to beat the noise_sd used below.
INPUT2D="$WORKDIR/input_2d"
mkdir -p "$INPUT2D"
python3 - "$INPUT2D" <<'EOF'
import json
import math
import random
import sys

out_dir = sys.argv[1]
N = 60
N_REF = 300
NOISE_SD = 0.05


def make_line_dataset(dataset_id, seed):
    random.seed(seed)
    coords = []
    for _ in range(N):
        t = random.uniform(-1.0, 1.0)
        coords.append([
            t + random.gauss(0.0, NOISE_SD),
            2.0 * t + 1.0 + random.gauss(0.0, NOISE_SD),
        ])
    ref_points = []
    for i in range(N_REF):
        t = -1.0 + 2.0 * i / (N_REF - 1)
        ref_points.append([t, 2.0 * t + 1.0])
    norm = math.sqrt(1.0 + 4.0)
    unit = [1.0 / norm, 2.0 / norm]
    ref_bases = [[unit] for _ in range(N_REF)]
    return {
        "dataset_id": dataset_id,
        "ambient_dim": 2,
        "generation_seed": seed,
        "manifolds": [{
            "manifold_id": 0,
            "intrinsic_dim": 1,
            "noise_sd": NOISE_SD,
            "reference_points": ref_points,
            "reference_tangent_bases": ref_bases,
        }],
        "points": {
            "coordinates": coords,
            "true_manifold_id": [0] * N,
        },
    }


for name, seed in [("line_a", 1), ("line_b", 2)]:
    with open(f"{out_dir}/{name}.json", "w") as fh:
        json.dump(make_line_dataset(name, seed), fh)
EOF

# A 3D counterpart, for the ambient_dim-mismatch rejection check below.
INPUT_MIXED="$WORKDIR/input_mixed"
mkdir -p "$INPUT_MIXED"
cp "$INPUT2D/line_a.json" "$INPUT_MIXED/"
python3 - "$INPUT_MIXED" <<'EOF'
import json
import random
import sys

out_dir = sys.argv[1]
random.seed(3)
N = 60
coords = [[random.gauss(0.0, 1.0), random.gauss(0.0, 1.0), random.gauss(0.0, 1.0)] for _ in range(N)]
dataset = {
    "dataset_id": "line_3d",
    "ambient_dim": 3,
    "manifolds": [{
        "manifold_id": 0,
        "intrinsic_dim": 1,
        "noise_sd": 0.05,
        "reference_points": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        "reference_tangent_bases": [[[0.577, 0.577, 0.577]], [[0.577, 0.577, 0.577]]],
    }],
    "points": {"coordinates": coords, "true_manifold_id": [0] * N},
}
with open(f"{out_dir}/line_3d.json", "w") as fh:
    json.dump(dataset, fh)
EOF

# An empty directory (no *.json files at all), for the "nothing to load" rejection check.
EMPTY_DIR="$WORKDIR/empty"
mkdir -p "$EMPTY_DIR"

run() {
  LD_LIBRARY_PATH="$ROOT/build" "$BIN" "$@"
}

failures=0
check() {
  if [[ "$1" -eq 0 ]]; then
    echo "ok: $2"
  else
    echo "FAIL: $2"
    failures=$((failures + 1))
  fi
}

echo "=== basic multi-dataset run ==="
OUT1="$WORKDIR/out1"
mkdir -p "$OUT1"
run --input-dir "$INPUT2D" --output-dir "$OUT1" --seed 42 --chains 1 --max-iterations 3 \
    --o 5 >/dev/null
check "$?" "basic run exits 0"
for f in line_a_chain00.csv line_a_chain01.csv line_b_chain00.csv line_b_chain01.csv; do
  test -s "$OUT1/$f"
  check "$?" "basic run wrote a non-empty $f"
done
head -1 "$OUT1/line_a_chain00.csv" | \
  grep -q "^iteration,k_min,chordal_dist_max_as_prcnt_of_range,d_max,G_max,RMSE_change_max,radius_percentile,min_stable_iterations,L_A,L_D,L_C,L_E,F,L,T,accepted$"
check "$?" "chain log has the expected CSV header"
LINES=$(wc -l < "$OUT1/line_a_chain00.csv")
[[ "$LINES" -eq 4 ]] # header + 3 iterations
check "$?" "chain log has one header row plus one row per iteration"

echo "=== --store-iteration-shatter-results ==="
OUT2="$WORKDIR/out2"
mkdir -p "$OUT2"
run --input-dir "$INPUT2D" --output-dir "$OUT2" --seed 7 --chains 0 --max-iterations 2 \
    --o 5 --store-iteration-shatter-results >/dev/null
check "$?" "--store-iteration-shatter-results run exits 0"
for iter in iter_000000 iter_000001; do
  DIR="$OUT2/line_a_chain00/$iter"
  test -d "$DIR"
  check "$?" "shatter output directory $DIR exists"
  for f in results.json points.csv ensemble_overlap_coefficients.csv super_ensembles.tsv; do
    test -f "$DIR/$f"
    check "$?" "shatter output $iter/$f was written"
  done
done
python3 -c "import json; json.load(open('$OUT2/line_a_chain00/iter_000000/results.json'))" \
  && shatter_json_ok=0 || shatter_json_ok=1
check "$shatter_json_ok" "shatter results.json is valid JSON"

echo "=== validation errors ==="
if run --output-dir "$OUT1" --seed 1 --o 5 >/dev/null 2>&1; then
  echo "FAIL: a run missing --input-dir should have been rejected"; failures=$((failures + 1))
else
  echo "ok: a run missing --input-dir rejected"
fi

if run --input-dir "$INPUT2D" --output-dir "$OUT1" --o 5 >/dev/null 2>&1; then
  echo "FAIL: a run missing --seed should have been rejected"; failures=$((failures + 1))
else
  echo "ok: a run missing --seed rejected"
fi

if run --input-dir "$INPUT2D" --output-dir "$OUT1" --seed 1 >/dev/null 2>&1; then
  echo "FAIL: a run missing --o should have been rejected"; failures=$((failures + 1))
else
  echo "ok: a run missing --o rejected"
fi

if run --input-dir "$EMPTY_DIR" --output-dir "$OUT1" --seed 1 --o 5 >/dev/null 2>&1; then
  echo "FAIL: an --input-dir with no *.json files should have been rejected"; failures=$((failures + 1))
else
  echo "ok: an --input-dir with no *.json files rejected"
fi

if run --input-dir "$INPUT_MIXED" --output-dir "$OUT1" --seed 1 --o 5 >/dev/null 2>&1; then
  echo "FAIL: an --input-dir mixing ambient_dim 2 and 3 should have been rejected"; failures=$((failures + 1))
else
  echo "ok: an --input-dir mixing ambient_dim 2 and 3 rejected"
fi

if [[ "$failures" -eq 0 ]]; then
  echo "ALL PASSED"
  exit 0
else
  echo "$failures check(s) FAILED"
  exit 1
fi
