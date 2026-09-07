#!/usr/bin/env python3
"""Generate independent reference events for the membrane/peptide fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
from MDAnalysis.lib.distances import distance_array, minimize_vectors


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation" / "data" / "membrane-peptide"
OUTPUT = ROOT / "validation" / "membrane_peptide_reference.json"
EXPECTED_VERSION = "2.10.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def intervals(times: list[float], states: list[bool]) -> list[list[float]]:
    result: list[list[float]] = []
    start: float | None = None
    previous = 0.0
    for time, state in zip(times, states, strict=True):
        if state and start is None:
            start = time
        elif not state and start is not None:
            result.append([start, previous])
            start = None
        previous = time
    if start is not None:
        result.append([start, previous])
    return result


def event_check(
    *,
    identifier: str,
    observable: str,
    implementation: str,
    query: str,
    times: list[float],
    states: list[bool],
    threshold: dict[str, object],
    clearance: float,
) -> dict[str, object]:
    found = intervals(times, states)
    return {
        "id": identifier,
        "observable": observable,
        "reference_implementation": implementation,
        "query": query,
        "threshold": threshold,
        "minimum_threshold_clearance": clearance,
        "matching_frame_count": sum(states),
        "event_count": len(found),
        "intervals_ps": found,
        "time_tolerance_ps": 0.001,
    }


def build_reference() -> dict[str, object]:
    import MDAnalysis as mda

    if mda.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            f"Reference generation requires MDAnalysis {EXPECTED_VERSION}; "
            f"received {mda.__version__}"
        )

    topology_path = DATA / "topology.pdb"
    trajectory_path = DATA / "trajectory.xtc"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Element information is missing.*")
        universe = mda.Universe(str(topology_path), str(trajectory_path))
    peptide = universe.select_atoms("chainID P")
    surface = universe.select_atoms("chainID M and name P")
    origin = universe.select_atoms("chainID P and resid 1 and name CA")
    target = universe.select_atoms("chainID P and resid 19 and name CA")
    if len(peptide) != 197 or len(surface) != 128:
        raise RuntimeError("Unexpected peptide or membrane-surface selection")
    if len(origin) != 1 or len(target) != 1:
        raise RuntimeError("Unexpected peptide endpoint selections")

    times: list[float] = []
    surface_distances: list[float] = []
    orientations: list[float] = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Reload offsets from trajectory.*")
        for timestep in universe.trajectory:
            if timestep.dimensions is None:
                raise RuntimeError("The membrane fixture must declare a periodic cell")
            times.append(round(float(timestep.time), 4))
            surface_distances.append(
                float(
                    distance_array(
                        peptide.positions, surface.positions, box=timestep.dimensions
                    ).min()
                    / 10.0
                )
            )
            vector = minimize_vectors(
                target.positions[0] - origin.positions[0], timestep.dimensions
            )
            cosine = abs(float(vector[2])) / float(np.linalg.norm(vector))
            orientations.append(
                float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))
            )

    surface_threshold = 0.281
    orientation_threshold = 30.0
    checks = [
        event_check(
            identifier="peptide_surface_distance_lt_0_281_nm",
            observable="SURFACE_DISTANCE",
            implementation=(
                "minimum MDAnalysis.lib.distances.distance_array value with the "
                "per-frame periodic cell"
            ),
            query="FIND SURFACE_DISTANCE(chain P, chain M and name P) < 0.281nm;",
            times=times,
            states=[value < surface_threshold for value in surface_distances],
            threshold={"comparison": "<", "value": surface_threshold, "unit": "nm"},
            clearance=float(
                np.min(np.abs(np.asarray(surface_distances) - surface_threshold))
            ),
        ),
        event_check(
            identifier="peptide_end_to_end_orientation_lt_30_deg",
            observable="ORIENTATION",
            implementation=(
                "MDAnalysis.lib.distances.minimize_vectors followed by a NumPy "
                "unoriented angle to Cartesian z"
            ),
            query=(
                "FIND ORIENTATION(chain P and resid 1 and name CA, "
                "chain P and resid 19 and name CA, axis=z) < 30deg;"
            ),
            times=times,
            states=[value < orientation_threshold for value in orientations],
            threshold={
                "comparison": "<",
                "value": orientation_threshold,
                "unit": "degree",
            },
            clearance=float(
                np.min(np.abs(np.asarray(orientations) - orientation_threshold))
            ),
        ),
    ]

    return {
        "schema_version": 1,
        "dataset": {
            "name": "peptide in a pure DMPC membrane",
            "source": "MDAnalysisData membrane_peptide",
            "source_url": "https://figshare.com/articles/dataset/Molecular_dynamics_trajectory_for_membrane_peptide_tutorial/8046437",
            "citation_doi": "10.6084/m9.figshare.8046437",
            "license": "CC BY 4.0",
            "source_files": {
                "source.tpr": "677a3ae55e35c24f37f2610eafa92d19285d1774731d6ffb9a99dfde39b8c437",
                "source.xtc": "f9bdfee4e1aa69ccfeef21cb74703202f6728f514543c4125382bd5250773eb7",
            },
            "derivation": "protein plus DMPC phosphorus atoms; every tenth source frame",
            "frame_count": len(times),
            "frame_spacing_ps": 100.0,
            "atom_count": len(universe.atoms),
        },
        "reference_software": {
            "MDAnalysis": EXPECTED_VERSION,
            "NumPy": np.__version__,
        },
        "prepared_files": {
            name: sha256(DATA / name) for name in ("topology.pdb", "trajectory.xtc")
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="fail if the committed reference is stale"
    )
    args = parser.parse_args()
    generated = build_reference()
    if args.check:
        recorded = json.loads(OUTPUT.read_text())
        if recorded != generated:
            raise SystemExit(
                "The recorded membrane/peptide reference is stale; regenerate it "
                "without --check"
            )
        print(f"Reference manifest is current: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(generated, indent=2) + "\n")
    print(f"Wrote independent reference manifest: {OUTPUT}")


if __name__ == "__main__":
    main()
