#!/usr/bin/env python3
"""
Combines interactive_template.html + vendor/d3.v7.min.js + <prefix>.json into one
self-contained, directly double-clickable <prefix>_interactive.html -- no local HTTP server
needed (plain file:// fetch() of a sibling JSON file is blocked by browsers' CORS rules for
local files, so the data is inlined at generation time instead of fetched at view time).

Run: python3 render_interactive.py <prefix>
  (re-)writes <prefix>.json via export_json.py if it does not already exist, then writes
  <prefix>_interactive.html.

interactive_template.html is the reusable part -- edit it to change the visualization: the
only thing that changes per report is the DATA JSON it gets built with.
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_PATH = os.path.join(HERE, "interactive_template.html")
D3_PATH = os.path.join(HERE, "vendor", "d3.v7.min.js")


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 render_interactive.py <prefix>", file=sys.stderr)
        sys.exit(1)
    prefix = sys.argv[1]

    json_path = f"{prefix}.json"
    if not os.path.exists(json_path):
        subprocess.run([sys.executable, os.path.join(HERE, "export_json.py"), prefix], check=True)

    with open(json_path) as fh:
        data = json.load(fh)
    with open(TEMPLATE_PATH) as fh:
        template = fh.read()
    with open(D3_PATH) as fh:
        d3_source = fh.read()

    # </script> inside the vendored D3 source or the JSON payload would prematurely close the
    # <script> tag it's embedded in -- escape the one substring that can trigger this.
    d3_safe = d3_source.replace("</script", "<\\/script")
    data_json = json.dumps(data).replace("</script", "<\\/script")

    html = template.replace("__STC_D3_JS__", d3_safe).replace("__STC_DATA__", data_json)

    out_path = f"{prefix}_interactive.html"
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"[render_interactive] wrote {out_path} ({len(html) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
