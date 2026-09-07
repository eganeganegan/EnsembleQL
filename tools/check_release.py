#!/usr/bin/env python3
"""Check that all public EnsembleQL versions agree with a release tag."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def python_version() -> str:
    tree = ast.parse((ROOT / "python" / "ensembleql" / "__init__.py").read_text())
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if (
            isinstance(target, ast.Name)
            and target.id == "__version__"
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            return node.value.value
    raise RuntimeError("Could not find ensembleql.__version__")


def cmake_version() -> str:
    text = (ROOT / "CMakeLists.txt").read_text()
    match = re.search(r"project\(EnsembleQL\s+VERSION\s+([^\s)]+)", text)
    if match is None:
        raise RuntimeError("Could not find the CMake project version")
    return match.group(1)


def citation_version() -> str:
    text = (ROOT / "CITATION.cff").read_text()
    match = re.search(r'^version:\s*["\']?([^"\'\n]+)', text, re.MULTILINE)
    if match is None:
        raise RuntimeError("Could not find the CITATION.cff version")
    return match.group(1).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tag", help="also require an exact release tag such as v0.1.0"
    )
    args = parser.parse_args()

    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    expected = project["project"]["version"]
    versions = {
        "pyproject.toml": expected,
        "python/ensembleql/__init__.py": python_version(),
        "CMakeLists.txt": cmake_version(),
        "CITATION.cff": citation_version(),
    }
    mismatches = {path: value for path, value in versions.items() if value != expected}
    if mismatches:
        details = ", ".join(f"{path}={value}" for path, value in mismatches.items())
        raise SystemExit(f"Release version mismatch; expected {expected}: {details}")
    if args.tag is not None and args.tag != f"v{expected}":
        raise SystemExit(
            f"Release tag {args.tag!r} does not match project version v{expected}"
        )
    changelog = (ROOT / "CHANGELOG.md").read_text()
    if f"## [{expected}]" not in changelog:
        raise SystemExit(f"CHANGELOG.md has no section for version {expected}")
    print(f"Release metadata agrees on version {expected}")


if __name__ == "__main__":
    main()
