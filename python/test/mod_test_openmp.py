#!/usr/bin/env python3
"""max_omp_threads() reports what the build does with OpenMP: 1 in a serial build, and in an
OpenMP build the runtime's maximum, which follows OMP_NUM_THREADS.

What the build was MEANT to be comes from TOX_EXPECT_OPENMP, which run_all_tests.sh sets from its
--openmp option (1 with it, 0 without). Run by hand without the variable, the test accepts either
answer, as it cannot tell which build it runs against."""

import os
import subprocess
import sys

PYTHON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, PYTHON_DIR)
from tensor_omics import max_omp_threads
from test_helpers import run_all_tests

N_THREADS_REQUESTED = 3

# TOX_EXPECT_OPENMP -> the answers max_omp_threads() may give under OMP_NUM_THREADS=3
ALLOWED_ANSWERS = {"1": (N_THREADS_REQUESTED,), "0": (1,), "": (1, N_THREADS_REQUESTED)}


def test_max_omp_threads_is_a_positive_int():
    max_threads = max_omp_threads()
    assert isinstance(max_threads, int), f"an int, not {type(max_threads).__name__}"
    assert max_threads >= 1, f"at least one thread, not {max_threads}"


def test_max_omp_threads_follows_omp_num_threads():
    expect_openmp = os.environ.get("TOX_EXPECT_OPENMP", "")
    assert expect_openmp in ALLOWED_ANSWERS, f"TOX_EXPECT_OPENMP is 1, 0 or unset, not {expect_openmp!r}"
    # in a new process: the OpenMP runtime reads OMP_NUM_THREADS once, when it is loaded, and
    # this one may have loaded it before this test could set the variable
    child = subprocess.run(
        [sys.executable, "-c", "from tensor_omics import max_omp_threads; print(max_omp_threads())"],
        env={**os.environ, "OMP_NUM_THREADS": str(N_THREADS_REQUESTED), "PYTHONPATH": PYTHON_DIR},
        capture_output=True, text=True, check=True,
    )
    max_threads = int(child.stdout)
    allowed = ALLOWED_ANSWERS[expect_openmp]
    assert max_threads in allowed, \
        f"{' or '.join(map(str, allowed))} with TOX_EXPECT_OPENMP={expect_openmp!r} and " \
        f"OMP_NUM_THREADS={N_THREADS_REQUESTED}, not {max_threads}"


if __name__ == "__main__":
    run_all_tests(globals().values())
