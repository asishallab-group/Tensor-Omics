#!/bin/bash
# Smoke test for stc_cli: builds the project, runs the CLI end to end (manual parameters,
# --estimate-parameters, and two validation-error paths), and checks the expected output
# files/exit codes. Not part of `test_runner.sh`'s Fortran/Python/R suites -- this exercises
# the C executable itself, which fpm builds but does not otherwise test.
#
# Usage: C-layer/test_stc_cli.sh [--skip-build]
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$PWD"

if [[ "${1:-}" != "--skip-build" ]]; then
  TOX_SKIP_CODE_GENERATION=1 ./build.sh
fi

BIN=$(find build -name "stc_cli" -type f -executable | head -1)
if [[ -z "$BIN" ]]; then
  echo "test_stc_cli.sh: could not find the built stc_cli executable" >&2
  exit 1
fi

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

FIXTURE="$WORKDIR/fixture.csv"
python3 - "$FIXTURE" <<'EOF'
import random
import sys
random.seed(42)
lines = ["x,y"]
for _ in range(15):
    lines.append(f"{random.gauss(0, 0.3):.4f},{random.gauss(0, 0.3):.4f}")
for _ in range(15):
    lines.append(f"{random.gauss(10, 0.3):.4f},{random.gauss(10, 0.3):.4f}")
with open(sys.argv[1], "w") as f:
    f.write("\n".join(lines) + "\n")
EOF

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

echo "=== manual parameters ==="
OUT1="$WORKDIR/out1"
mkdir -p "$OUT1"
run --input "$FIXTURE" --header --n-records 30 --output-dir "$OUT1" \
    --k-min 5 --k-density 5 --chordal-dist-max-as-prcnt-of-range 0.3 --d-max 1 --g-max 5.0 \
    --rmse-change-max 5.0 --o 4 --exclusion-radius-percentile 60 >/dev/null
check "$?" "manual run exits 0"
for f in report.html results.json points.csv ensemble_overlap_coefficients.csv super_ensembles.tsv; do
  test -s "$OUT1/$f"
  check "$?" "manual run wrote a non-empty $f"
done
grep -q '"dim_names":\["x","y"\]' "$OUT1/results.json"
check "$?" "manual run's JSON captured the CSV header as dim_names"
head -1 "$OUT1/points.csv" | grep -q "^row,ensembles,super_ensembles,low_confidence_ensembles,seed_of$"
check "$?" "manual run's points.csv has the expected header"

echo "=== --estimate-parameters ==="
OUT2="$WORKDIR/out2"
mkdir -p "$OUT2"
run --input "$FIXTURE" --header --n-records 30 --output-dir "$OUT2" \
    --estimate-parameters --rmse-change-max 5.0 --o 4 --seed-max-set-size 30 --n-anchors 5 >/dev/null
check "$?" "--estimate-parameters run exits 0"
grep -q '"estimated_k_min"' "$OUT2/results.json"
check "$?" "--estimate-parameters run reports estimated_k_min in the JSON"

echo "=== validation errors ==="
if run --input "$FIXTURE" --header --n-records 30 --output-dir "$OUT1" \
       --estimate-parameters --k-min 5 --rmse-change-max 5.0 --o 4 >/dev/null 2>&1; then
  echo "FAIL: --estimate-parameters + --k-min should have been rejected"
  failures=$((failures + 1))
else
  echo "ok: --estimate-parameters + --k-min rejected"
fi

if run --input "$FIXTURE" --header --n-records 30 --output-dir "$OUT1" >/dev/null 2>&1; then
  echo "FAIL: a run missing required flags should have been rejected"
  failures=$((failures + 1))
else
  echo "ok: a run missing required flags rejected"
fi

if [[ "$failures" -eq 0 ]]; then
  echo "ALL PASSED"
  exit 0
else
  echo "$failures check(s) FAILED"
  exit 1
fi
