#!/usr/bin/env python3
"""Prepare the pinned MDAnalysis adenylate-kinase validation fixture."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from shutil import copy2
import warnings


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation" / "data" / "adk"
EXPECTED_VERSION = "2.10.0"
EXPECTED_SOURCE_SHA256 = {
    "adk.psf": "96cec916c4b5b19a7acb91bd1e672fb9c8b032aee1e7f0868a50655902c9b5b8",
    "adk_dims.dcd": "859a5bd9e7de45a0f2401f7971c5ffc296168e7c23f6f65382bde7c1686c19f1",
}
EXPECTED_OUTPUT_SHA256 = {
    "topology.pdb": "7bca72a84d8606d1091e2d9cb632563210db16d00afe73ec77f8c73dcb945ea6",
    "trajectory.dcd": "859a5bd9e7de45a0f2401f7971c5ffc296168e7c23f6f65382bde7c1686c19f1",
    "trajectory.xtc": "aa7b9a6e70b5b13719d60193ad5773330e4df851f173d4129a4de9a46b11ad4f",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_digest(path: Path, expected: str) -> None:
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(
            f"SHA-256 mismatch for {path}: expected {expected}, received {actual}"
        )


def check_outputs() -> None:
    for name, expected in EXPECTED_OUTPUT_SHA256.items():
        path = OUTPUT / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing prepared validation file: {path}")
        require_digest(path, expected)


def prepare(force: bool) -> None:
    import MDAnalysis as mda
    import MDAnalysisTests
    from MDAnalysisTests.datafiles import DCD, PSF

    if mda.__version__ != EXPECTED_VERSION or MDAnalysisTests.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            "Validation preparation requires MDAnalysis and MDAnalysisTests "
            f"{EXPECTED_VERSION}; received {mda.__version__} and "
            f"{MDAnalysisTests.__version__}"
        )

    sources = {"adk.psf": Path(PSF), "adk_dims.dcd": Path(DCD)}
    for name, path in sources.items():
        require_digest(path, EXPECTED_SOURCE_SHA256[name])

    if not force and all((OUTPUT / name).is_file() for name in EXPECTED_OUTPUT_SHA256):
        check_outputs()
        print(f"Validation fixture is already prepared in {OUTPUT}")
        return

    OUTPUT.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe(PSF, DCD)
    universe.trajectory[0]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unit cell dimensions not found.*")
        warnings.filterwarnings("ignore", message="Found no information for attr:.*")
        warnings.filterwarnings("ignore", message="Found missing chainIDs.*")
        universe.atoms.write(OUTPUT / "topology.pdb", bonds="all")
        universe.atoms.write(OUTPUT / "trajectory.xtc", frames="all")
    copy2(DCD, OUTPUT / "trajectory.dcd")

    check_outputs()
    print(f"Prepared validated 98-frame fixture in {OUTPUT}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="verify existing outputs without writing"
    )
    parser.add_argument(
        "--force", action="store_true", help="replace existing generated outputs"
    )
    args = parser.parse_args()
    if args.check:
        check_outputs()
        print(f"Validation fixture checksums match in {OUTPUT}")
    else:
        prepare(args.force)


if __name__ == "__main__":
    main()
