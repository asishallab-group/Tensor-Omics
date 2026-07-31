"""Command-line entry point.

A thin shell over `generate`: parse arguments, run the pipeline, print diagnostics, set an
exit code. All the work is in `generate.py`, so this stays small enough to read at a glance.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import Paths
from .generate import generate, generate_and_write

C_BINDING_TARGETS = ("c", "python", "r")
SNIPPETS_TARGETS = ("snippets",)
ALL_TARGETS = (*C_BINDING_TARGETS, *SNIPPETS_TARGETS)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="generate_code",
        description="Generate the C, Python and R bindings from the Fortran sources.",
    )
    parser.add_argument(
        "--root", type=Path, default=Path("."),
        help="the repository root; defaults to the current directory",
    )
    parser.add_argument(
        "--src", type=Path, default=None,
        help="the source directory to read, relative to root (default: src)",
    )
    parser.add_argument(
        "--target", choices=ALL_TARGETS, action="append", dest="targets",
        help=f"a target to generate; repeatable, defaults to ({", ".join(C_BINDING_TARGETS)})",
    )
    parser.add_argument(
        "--library", default="build/libtensor-omics.so",
        help="where the built shared library is, as the Python loader should look for it",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="report problems without writing anything, and exit non-zero if any files "
             "would change",
    )
    parser.add_argument(
        "--no-clean", action="store_true",
        help="do not remove the output directories first (may leave stale files behind)",
    )
    parser.add_argument(
        "--no-color", action="store_true", help="do not colour the diagnostics",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    paths = Paths(root=args.root)
    if args.src is not None:
        paths = Paths(root=args.root, src_dir=args.src)
    targets = tuple(args.targets) if args.targets else C_BINDING_TARGETS
    color = not args.no_color and sys.stderr.isatty()

    if args.check:
        result = generate(paths, targets, library=args.library)
    else:
        result = generate_and_write(
            paths, targets, library=args.library, clean=not args.no_clean
        )

    diagnostics = result.diagnostics
    if len(diagnostics):
        print(diagnostics.render(color=color), file=sys.stderr)

    if not result.ok:
        print(
            f"\n{len(diagnostics.errors)} error(s); nothing was written.", file=sys.stderr
        )
        return 1

    if args.check:
        changed = _changed(result)
        if changed:
            print(f"{len(changed)} file(s) would change:", file=sys.stderr)
            for path in changed:
                print(f"  {path}", file=sys.stderr)
            return 1
        print("up to date.", file=sys.stderr)
        return 0

    print(
        f"generated {len(result.files)} file(s) "
        f"({len(diagnostics.warnings)} warning(s)).",
        file=sys.stderr,
    )
    return 0


def _changed(result) -> list[Path]:
    """Which generated files differ from what is on disk, for --check."""
    changed = []
    for file in result.files:
        try:
            existing = file.path.read_text()
        except OSError:
            existing = None
        if existing != file.content:
            changed.append(file.path)
    return changed


if __name__ == "__main__":
    raise SystemExit(main())
