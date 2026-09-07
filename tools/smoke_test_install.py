#!/usr/bin/env python3
"""Smoke-test an installed release wheel outside the source package path."""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

import ensembleql as eql


def main() -> None:
    if eql.__version__ != version("ensembleql"):
        raise RuntimeError("Package and distribution versions disagree")
    if not eql.chemfiles_backend_available():
        raise RuntimeError("Release wheels must include the Chemfiles backend")

    root = Path(__file__).resolve().parents[1]
    example = root / "examples" / "idr_contact_switching"
    trajectory = eql.load(
        example / "switching.xyz", topology=example / "switching.pdb"
    )
    events = trajectory.query("FIND CONTACT(resid 17, resid 42);")
    if len(events) != 1 or events[0].start != 0.0 or events[0].end != 2000.0:
        raise RuntimeError("Installed-wheel query smoke test returned unexpected events")
    print(f"Validated installed EnsembleQL {eql.__version__} release wheel")


if __name__ == "__main__":
    main()
