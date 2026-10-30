#!/usr/bin/env python3
"""max_omp_threads() reports what the build does with OpenMP: 1 in a serial build, and in an
OpenMP build the runtime's maximum, which follows OMP_NUM_THREADS."""

import os
import subprocess
import sys

PYTHON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, PYTHON_DIR)
from tensor_omics import max_omp_threads
from test_helpers import run_all_tests


def test_max_omp_threads_is_a_positive_int():
    max_threads = max_omp_threads()
    assert isinstance(max_threads, int), f"an int, not {type(max_threads).__name__}"
    assert max_threads >= 1, f"at least one thread, not {max_threads}"


def test_max_omp_threads_follows_omp_num_threads():
    # in a new process: the OpenMP runtime reads OMP_NUM_THREADS once, when it is loaded, and
    # this one may have loaded it before this test could set the variable
    child = subprocess.run(
        [sys.executable, "-c", "from tensor_omics import max_omp_threads; print(max_omp_threads())"],
        env={**os.environ, "OMP_NUM_THREADS": "3", "PYTHONPATH": PYTHON_DIR},
        capture_output=True, text=True, check=True,
    )
    max_threads = int(child.stdout)
    assert max_threads in (1, 3), f"1 (serial) or OMP_NUM_THREADS=3 (OpenMP), not {max_threads}"


if __name__ == "__main__":
    run_all_tests(globals().values())
