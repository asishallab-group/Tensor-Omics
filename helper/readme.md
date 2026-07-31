# helper

Developer tooling for TensorOmics. Not compiled into the library.

| Path | What it is |
|---|---|
| [`codegen/`](codegen/) | the interface generator: turns exported Fortran procedures into the C, Python and R interfaces. See its [README](codegen/README.md). |
| `generate_code.py` | entry point for the generator (`python helper/generate_code.py --help`); also emits the VS Code snippets (`--target snippets`) |
