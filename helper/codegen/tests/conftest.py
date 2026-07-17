"""Test configuration.

Puts `helper/` on the import path so tests import the generator the same way
`helper/generate_code.py` does, i.e. as the top-level package `codegen`.
"""

import sys
from pathlib import Path

HELPER_DIR = Path(__file__).resolve().parents[2]
REPO_ROOT = HELPER_DIR.parent

if str(HELPER_DIR) not in sys.path:
    sys.path.insert(0, str(HELPER_DIR))
