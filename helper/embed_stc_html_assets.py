#!/usr/bin/env python3
"""
One-off generator (NOT part of the main codegen pipeline in helper/codegen/): reads the
vendored D3 and Babylon.js bundles and the interactive report template from
misc/STC-experiments/, and emits src/tox_stc_html_assets.F90, a generated (never hand-edited)
Fortran module holding their content as compile-time character parameters -- so the resulting
Fortran/C/Python/R binary/package never needs to locate these files at runtime; they are baked
in. D3 backs the 2D report, Babylon.js backs the 3+D report (see misc/mod_STC.md,
"Visualization for 3+ dimensions") -- both ship in every build since which one a given report
needs is only known at report-generation time (from the run's own n_dimensions), not at
compile time.

Run: python3 helper/embed_stc_html_assets.py
Re-run this whenever misc/STC-experiments/vendor/d3.v7.min.js,
misc/STC-experiments/vendor/babylon.js, or misc/STC-experiments/interactive_template.html
change (e.g. a version bump) -- the generated .F90 file is committed, like everything under
src/generated/.

Why compile-time embedding, not a runtime-read external file: this content ships inside
Tensor Omics regardless of delivery shape (a bare Fortran+C executable, a Python wheel, an R
package) -- runtime path resolution (relative to argv[0], an installed package's own data
directory, an environment variable, ...) is a different, fragile answer per delivery
mechanism. Baking the content into the compiled artifact itself sidesteps that class of
problem entirely: there is nothing to "find" at runtime.

Why a hand-rolled generator instead of e.g. a C `xxd -i`-style .inc: this needs to become a
valid Fortran CHARACTER parameter, quote-escaped, chunked to stay under this build's
practical free-form line-length limit (empirically well under 16000 chars/line -- one
~4000-char chunk per generated source line stays comfortably inside that; see the commit
introducing this file for the original, much smaller-chunk experiment). Chunk size matters
well beyond the line-length ceiling: gfortran's own handling of a `//`-concatenation
expression chain scales badly with the number of pieces, not their total length -- the
original 90-char/3-per-line grouping (fine for D3's ~280KB) took several minutes to compile
once this file started also embedding Babylon.js's ~8.3MB single-line minified bundle at the
same granularity (tens of thousands of pieces); widening to ~4000 chars/piece is a ~44x
reduction in piece count for the same content, with no observed downside. Newlines are
represented via achar(10) rather than embedded raw in the generated source (a raw newline
byte inside a Fortran string literal would terminate the *generated source's own* line,
corrupting it).
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)

D3_PATH = os.path.join(REPO_ROOT, "misc", "STC-experiments", "vendor", "d3.v7.min.js")
BABYLON_PATH = os.path.join(REPO_ROOT, "misc", "STC-experiments", "vendor", "babylon.js")
TEMPLATE_PATH = os.path.join(REPO_ROOT, "misc", "STC-experiments", "interactive_template.html")
OUT_PATH = os.path.join(REPO_ROOT, "src", "tox_stc_html_assets.F90")

D3_PLACEHOLDER = "__STC_D3_JS__"
BABYLON_PLACEHOLDER = "__STC_BABYLON_JS__"
DATA_PLACEHOLDER = "__STC_DATA__"

RAW_CHUNK_SIZE = 4000
CHUNKS_PER_LINE = 1
INDENT = " " * 8


def fortran_literal_pieces(text):
    """Split `text` into a flat sequence of Fortran expression pieces (each either a quoted
    string literal or the expression achar(10)) that, concatenated with `//`, reproduce
    `text` exactly. Newlines are extracted first and represented as achar(10) so no chunk
    ever contains a raw newline byte (which would break the *generated* source line)."""
    pieces = []
    lines = text.split("\n")
    for i, line in enumerate(lines):
        for start in range(0, len(line), RAW_CHUNK_SIZE):
            raw_chunk = line[start:start + RAW_CHUNK_SIZE]
            escaped = raw_chunk.replace("'", "''")
            pieces.append("'" + escaped + "'")
        if i < len(lines) - 1:
            pieces.append("achar(10)")
    if not pieces:
        pieces.append("''")
    return pieces


def render_parameter(name, text):
    pieces = fortran_literal_pieces(text)
    lines = [f"    character(len=*), parameter :: {name} = &"]
    for i in range(0, len(pieces), CHUNKS_PER_LINE):
        group = pieces[i:i + CHUNKS_PER_LINE]
        is_last = (i + CHUNKS_PER_LINE >= len(pieces))
        joined = "//".join(group)
        suffix = "" if is_last else "//&"
        lines.append(INDENT + joined + suffix)
    return "\n".join(lines)


def main():
    with open(D3_PATH, "r", encoding="utf-8") as fh:
        d3_js = fh.read()
    with open(BABYLON_PATH, "r", encoding="utf-8") as fh:
        babylon_js = fh.read()
    with open(TEMPLATE_PATH, "r", encoding="utf-8") as fh:
        template = fh.read()

    for placeholder in (D3_PLACEHOLDER, BABYLON_PLACEHOLDER, DATA_PLACEHOLDER):
        if placeholder not in template:
            raise ValueError(f"{TEMPLATE_PATH} is missing the {placeholder} placeholder")

    head, rest = template.split(D3_PLACEHOLDER, 1)
    mid1, rest = rest.split(BABYLON_PLACEHOLDER, 1)
    mid2, tail = rest.split(DATA_PLACEHOLDER, 1)

    parts = [
        render_parameter("REPORT_TEMPLATE_HEAD", head),
        render_parameter("D3_JS", d3_js),
        render_parameter("REPORT_TEMPLATE_MID1", mid1),
        render_parameter("BABYLON_JS", babylon_js),
        render_parameter("REPORT_TEMPLATE_MID2", mid2),
        render_parameter("REPORT_TEMPLATE_TAIL", tail),
    ]

    out = (
        "!> Generated by helper/embed_stc_html_assets.py from\n"
        "!| misc/STC-experiments/vendor/d3.v7.min.js,\n"
        "!| misc/STC-experiments/vendor/babylon.js, and\n"
        "!| misc/STC-experiments/interactive_template.html. Do not edit -- regenerate instead.\n"
        "!|\n"
        "!| The interactive HTML report is assembled by writing, in order:\n"
        "!| REPORT_TEMPLATE_HEAD, D3_JS, REPORT_TEMPLATE_MID1, BABYLON_JS, REPORT_TEMPLATE_MID2,\n"
        "!| <the run's JSON payload>, REPORT_TEMPLATE_TAIL -- see\n"
        "!| tox_stc_json::write_stc_interactive_html_report. D3 backs the 2D report, Babylon.js\n"
        "!| backs the 3+D report; both are always baked in since which one a given run needs is\n"
        "!| only known at report-generation time (n_dimensions), not at compile time.\n"
        "module tox_stc_html_assets\n"
        "    implicit none\n"
        "\n"
        "    private\n"
        "    public :: REPORT_TEMPLATE_HEAD, REPORT_TEMPLATE_MID1, REPORT_TEMPLATE_MID2, &\n"
        "        REPORT_TEMPLATE_TAIL, D3_JS, BABYLON_JS\n"
        "\n"
        + "\n\n".join(parts) + "\n"
        "\n"
        "end module tox_stc_html_assets\n"
    )

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        fh.write(out)
    print(f"[embed_stc_html_assets] wrote {OUT_PATH} ({len(out) / 1024:.0f} KB source, "
          f"D3 {len(d3_js)} bytes, Babylon {len(babylon_js)} bytes, template {len(template)} bytes)")


if __name__ == "__main__":
    main()
