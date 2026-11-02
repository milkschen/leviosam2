# Development and formatting

See [INSTALL.md](INSTALL.md) for the build and test environments.

## Reproducible formatting

Use Python 3.9 or newer to create a separate formatting environment. The pinned
Black and clang-format versions are shared by local commands and CI:

```sh
python3 -m venv .venv-format
. .venv-format/bin/activate
python -m pip install -r requirements-format.txt
python scripts/format_code.py          # apply formatting
python scripts/format_code.py --check  # CI check, no writes
```

The C++ configuration follows the existing Google-based layout with four-space
indentation and an 80-column limit. Include order and documentation wrapping are
preserved. Pointer/reference alignment is inferred from each file's existing
majority style (with left alignment as the fallback), avoiding a separate
pointer-style migration. Formatting is deterministic for the pinned tool and
current file contents. Black uses its standard 88-column Python layout.

The runner selects tracked first-party files using Git: C++ files under `src/`
and Python files throughout the repository. Imported headers/libraries, C
sources, generated/build files and untracked files are excluded. The explicit
upstream-file exclusions live in `scripts/format_code.py`; update that list when
adding imported code. Stage new first-party files before running the formatter
so they are included. Do not reformat vendored libraries when updating them.

Formatting checks run on every PR, including PRs adding new files. Keep
formatting-only changes separate from behavior changes. When updating formatter
versions, change the requirements and runner version checks, review resulting
diffs, and verify that a second formatting pass produces no further changes.
