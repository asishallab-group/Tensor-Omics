#!/usr/bin/env python3
"""Generate the C, Python and R interfaces from the Fortran sources.

Run from the repository root:

    python helper/generate_code.py            # generate everything
    python helper/generate_code.py --check    # report problems, write nothing
    python helper/generate_code.py --help     # all options

Everything lives in the `codegen` package next to this file; this is only the entry point.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from codegen.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
