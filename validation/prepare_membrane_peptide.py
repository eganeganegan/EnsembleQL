#!/usr/bin/env python3
"""Prepare a compact published membrane/peptide validation fixture."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen
import warnings

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "validation" / "data" / "membrane-peptide"
EXPECTED_VERSION = "2.10.0"
SOURCE = {
    "source.tpr": {
        "url": "https://ndownloader.figshare.com/files/14993171",
        "sha256": "677a3ae55e35c24f37f2610eafa92d19285d1774731d6ffb9a99dfde39b8c437",
    },
    "source.xtc": {
        "url": "https://ndownloader.figshare.com/files/14993174",
        "sha256": "f9bdfee4e1aa69ccfeef21cb74703202f6728f514543c4125382bd5250773eb7",
    },
}
EXPECTED_OUTPUT_SHA256 = {
    "topology.pdb": "6ce495ffbbf3db724c0b58926a32105160271c68feb672c4259864a50af05775",
    "trajectory.xtc": "b3f5482f8029abbd9f85ddf33ff48e1e6c72cd2bbc4460227fe120beb056ee95",
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


def download_sources() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, metadata in SOURCE.items():
        destination = OUTPUT / name
        if destination.is_file() and sha256(destination) == metadata["sha256"]:
            continue
        partial = destination.with_suffix(destination.suffix + ".part")
        with urlopen(metadata["url"]) as response, partial.open("wb") as stream:
            while chunk := response.read(1024 * 1024):
                stream.write(chunk)
        require_digest(partial, metadata["sha256"])
        partial.replace(destination)


def check_outputs() -> None:
    for name, expected in EXPECTED_OUTPUT_SHA256.items():
        path = OUTPUT / name
        if not path.is_file():
            raise FileNotFoundError(f"Missing prepared validation file: {path}")
        require_digest(path, expected)


def prepare(force: bool, download: bool) -> None:
    import MDAnalysis as mda

    if mda.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            f"Validation preparation requires MDAnalysis {EXPECTED_VERSION}; "
            f"received {mda.__version__}"
        )
    if download:
        download_sources()
    for name, metadata in SOURCE.items():
        path = OUTPUT / name
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing source file: {path}; rerun with --download"
            )
        require_digest(path, metadata["sha256"])

    if not force and all((OUTPUT / name).is_file() for name in EXPECTED_OUTPUT_SHA256):
        check_outputs()
        print(f"Validation fixture is already prepared in {OUTPUT}")
        return

    universe = mda.Universe(
        str(OUTPUT / "source.tpr"), str(OUTPUT / "source.xtc")
    )
    universe.atoms.chainIDs = np.where(
        universe.atoms.segids == "seg_0_Protein",
        "P",
        np.where(universe.atoms.segids == "seg_1_DMPC", "M", "W"),
    )
    selected = universe.select_atoms("protein or (resname DMPC and name P)")
    if len(selected) != 325:
        raise RuntimeError(f"Expected 325 selected atoms, received {len(selected)}")

    universe.trajectory[0]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Found no information for attr:.*")
        warnings.filterwarnings("ignore", message="Found missing chainIDs.*")
        selected.write(str(OUTPUT / "topology.pdb"), bonds="all")
    with mda.Writer(str(OUTPUT / "trajectory.xtc"), n_atoms=len(selected)) as writer:
        for _ in universe.trajectory[::10]:
            writer.write(selected)

    check_outputs()
    print(f"Prepared validated 101-frame fixture in {OUTPUT}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="verify existing outputs without writing"
    )
    parser.add_argument(
        "--download", action="store_true", help="download the pinned Figshare sources"
    )
    parser.add_argument(
        "--force", action="store_true", help="replace existing generated outputs"
    )
    args = parser.parse_args()
    if args.check:
        check_outputs()
        print(f"Validation fixture checksums match in {OUTPUT}")
    else:
        prepare(args.force, args.download)


if __name__ == "__main__":
    main()
