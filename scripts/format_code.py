#!/usr/bin/env python3
"""Format tracked first-party C++ and Python; --check verifies without writing."""

import argparse
from importlib.metadata import PackageNotFoundError, version
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
# Imported libraries retain their upstream formatting. Also exclude local
# modifications to imported code: formatting them would obscure upstream diffs.
VENDORED = {
    "src/IITree.h",
    "src/bam.h",
    "src/gzstream.C",
    "src/gzstream.h",
    "src/ksw2.h",
    "src/rapidyaml.hpp",
    "src/robin_hood.h",
}
EXCLUDED_DIRS = {".git", ".pixi", "build", "vendor", "third_party", "sse2neon"}


def source_files():
    """Select tracked sources only, excluding imported and generated code."""
    result = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, check=True, capture_output=True
    )
    cpp, python = [], []
    for name in os.fsdecode(result.stdout).split("\0"):
        path = Path(name)
        if not name or name in VENDORED or EXCLUDED_DIRS.intersection(path.parts):
            continue
        if not (ROOT / path).is_file():
            continue
        if path.suffix == ".py":
            python.append(name)
        elif path.parts[0] == "src" and path.suffix in {".cpp", ".hpp", ".h", ".C"}:
            cpp.append(name)
    return sorted(cpp), sorted(python)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check without rewriting")
    args = parser.parse_args()
    for package, expected in [("black", "24.10.0"), ("clang-format", "18.1.8")]:
        try:
            installed = version(package)
        except PackageNotFoundError:
            installed = "not installed"
        if installed != expected:
            parser.error(
                f"{package} {expected} required (found {installed}); "
                "install requirements-format.txt in this Python environment"
            )
    clang = Path(sys.executable).parent / "clang-format"
    if not clang.is_file():
        found = shutil.which("clang-format")
        if not found:
            parser.error("clang-format executable not found")
        clang = Path(found)
    actual = subprocess.check_output([str(clang), "--version"], text=True)
    if "version 18.1.8" not in actual:
        parser.error(f"Wrong clang-format executable: {actual.strip()}")
    cpp, python = source_files()
    failed = False
    if cpp:
        options = ["--dry-run", "--Werror"] if args.check else ["-i"]
        result = subprocess.run([str(clang), "--style=file", *options, *cpp], cwd=ROOT)
        failed |= result.returncode != 0
    if python:
        options = ["--check"] if args.check else []
        result = subprocess.run(
            [sys.executable, "-m", "black", "--workers", "1", *options, *python],
            cwd=ROOT,
        )
        failed |= result.returncode != 0
    print(
        f"{'Checked' if args.check else 'Formatted'} {len(cpp)} C++ and {len(python)} Python files."
    )
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
