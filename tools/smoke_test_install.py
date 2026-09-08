#!/usr/bin/env python3
"""Smoke-test an installed release wheel outside the source package path."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from importlib.metadata import version
import json
from pathlib import Path

import ensembleql as eql
from ensembleql.cli import main as cli_main


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

    transition = trajectory.query(
        "FIND CONTACT(resid 17, resid 42) "
        "FOLLOWED_BY CONTACT(resid 17, resid 53) WITHIN 3ns;"
    )
    if len(transition) != 1 or transition[0].end != 4000.0:
        raise RuntimeError("Installed-wheel temporal query returned unexpected events")

    output = StringIO()
    with redirect_stdout(output):
        status = cli_main(
            [
                "query",
                "--topology",
                str(example / "switching.pdb"),
                "--trajectory",
                str(example / "switching.xyz"),
                "--query",
                "FIND CONTACT(resid 17, resid 42);",
                "--format",
                "json",
            ]
        )
    records = json.loads(output.getvalue())
    if status != 0 or len(records) != 1 or records[0]["end_ps"] != 2000.0:
        raise RuntimeError("Installed-wheel CLI smoke test returned unexpected output")
    print(f"Validated installed EnsembleQL {eql.__version__} release wheel")


if __name__ == "__main__":
    main()
