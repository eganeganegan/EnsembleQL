#!/usr/bin/env python3
"""Generate independent reference events for the prepared AdK fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.lib.distances import calc_dihedrals, distance_array


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation" / "data" / "adk"
OUTPUT = ROOT / "validation" / "adk_reference.json"
EXPECTED_VERSION = "2.10.0"
MDTRAJ_VERSION = "1.11.1.post2"


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
    **metadata: object,
) -> dict[str, object]:
    found = intervals(times, states)
    return {
        "id": identifier,
        "observable": observable,
        "reference_implementation": implementation,
        "query": query,
        **metadata,
        "matching_frame_count": sum(states),
        "event_count": len(found),
        "intervals_ps": found,
        "time_tolerance_ps": 0.0001,
    }


def build_reference() -> dict[str, object]:
    import MDAnalysis as mda
    import MDAnalysisTests
    import mdtraj
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
    from MDAnalysisTests.datafiles import DCD, PSF

    if mda.__version__ != EXPECTED_VERSION or MDAnalysisTests.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            "Reference generation requires MDAnalysis and MDAnalysisTests "
            f"{EXPECTED_VERSION}"
        )
    if mdtraj.__version__ != MDTRAJ_VERSION:
        raise RuntimeError(
            f"Reference generation requires MDTraj {MDTRAJ_VERSION}; "
            f"received {mdtraj.__version__}"
        )

    topology_path = DATA / "topology.pdb"
    trajectory_path = DATA / "trajectory.xtc"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Element information is missing.*")
        universe = mda.Universe(str(topology_path), str(trajectory_path))

    distance_first = universe.select_atoms("resid 1 and name CA")
    distance_second = universe.select_atoms("resid 2 and name CA")
    contacts_first = universe.select_atoms("resid 30:59")
    contacts_second = universe.select_atoms("resid 120:159")
    protein = universe.select_atoms("protein")
    protein_ca = universe.select_atoms("protein and name CA")
    dihedral_atoms = [
        universe.select_atoms(selection)
        for selection in (
            "resid 49 and name C",
            "resid 50 and name N",
            "resid 50 and name CA",
            "resid 50 and name C",
        )
    ]
    ring_a = universe.select_atoms("resid 181 and (name CG or name CD1 or name CD2)")
    ring_b = universe.select_atoms("resid 182 and (name CG or name CD1 or name CD2)")
    if len(distance_first) != 1 or len(distance_second) != 1:
        raise RuntimeError("Expected exactly one C-alpha atom in residues 1 and 2")
    if len(protein) != len(universe.atoms) or len(protein_ca) != 214:
        raise RuntimeError("Unexpected protein selections in the AdK topology")
    if any(len(group) != 1 for group in dihedral_atoms):
        raise RuntimeError("Unexpected residue-50 dihedral selections")
    if len(ring_a) != 3 or len(ring_b) != 3:
        raise RuntimeError("Unexpected aromatic ring selections")

    backbone: list[tuple[object, object, object, object, object]] = []
    for resid in range(2, 214):
        groups = tuple(
            universe.select_atoms(selection)
            for selection in (
                f"resid {resid - 1} and name C",
                f"resid {resid} and name N",
                f"resid {resid} and name CA",
                f"resid {resid} and name C",
                f"resid {resid + 1} and name N",
            )
        )
        if all(len(group) == 1 for group in groups):
            backbone.append(groups)

    salt_positive = universe.select_atoms("resid 78 and name NH2")
    salt_negative = universe.select_atoms("resid 79 and name OD1")
    times: list[float] = []
    distances: list[float] = []
    minimum_contact_distances: list[float] = []
    atom_contact_counts: list[int] = []
    residue_contact_counts: list[int] = []
    coordination_numbers: list[float] = []
    unweighted_rg: list[float] = []
    mass_weighted_rg: list[float] = []
    fitted_rmsds: list[float] = []
    raw_rmsds: list[float] = []
    dihedrals: list[float] = []
    helix_fractions: list[float] = []
    salt_bridge_distances: list[float] = []
    aromatic_distances: list[float] = []
    aromatic_angles: list[float] = []
    reference_ca: np.ndarray | None = None
    reference_protein: np.ndarray | None = None

    for timestep in universe.trajectory:
        times.append(round(float(timestep.time), 4))
        distances.append(float(np.linalg.norm(distance_first.positions[0] - distance_second.positions[0]) / 10.0))

        pair_distances = distance_array(contacts_first.positions, contacts_second.positions)
        minimum_contact_distances.append(float(pair_distances.min() / 10.0))
        contacting = np.argwhere(pair_distances <= 4.5)
        atom_contact_counts.append(int(len(contacting)))
        residue_pairs = {
            (int(contacts_first[int(first)].resid), int(contacts_second[int(second)].resid))
            for first, second in contacting
        }
        residue_contact_counts.append(len(residue_pairs))
        coordination_numbers.append(float(len(contacting) / len(contacts_first)))

        coordinates = protein.positions.astype(np.float64)
        centered = coordinates - coordinates.mean(axis=0)
        unweighted_rg.append(float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))) / 10.0))
        masses = protein.masses.astype(np.float64)
        mass_center = np.average(coordinates, axis=0, weights=masses)
        mass_weighted_rg.append(
            float(np.sqrt(np.average(np.sum((coordinates - mass_center) ** 2, axis=1), weights=masses)) / 10.0)
        )

        if reference_ca is None:
            reference_ca = protein_ca.positions.copy()
            reference_protein = protein.positions.copy()
        fitted_rmsds.append(float(rmsd(protein_ca.positions, reference_ca, center=True, superposition=True) / 10.0))
        raw_rmsds.append(float(rmsd(protein.positions, reference_protein, center=False, superposition=False) / 10.0))

        dihedrals.append(
            float(np.degrees(calc_dihedrals(*(group.positions for group in dihedral_atoms))[0]))
        )
        helical = 0
        for previous_c, nitrogen, alpha_c, carbon, next_n in backbone:
            phi = np.degrees(calc_dihedrals(previous_c.positions, nitrogen.positions, alpha_c.positions, carbon.positions)[0])
            psi = np.degrees(calc_dihedrals(nitrogen.positions, alpha_c.positions, carbon.positions, next_n.positions)[0])
            helical += -100.0 <= phi <= -30.0 and -80.0 <= psi <= -5.0
        helix_fractions.append(float(helical / len(backbone)))

        salt_bridge_distances.append(float(np.linalg.norm(salt_positive.positions[0] - salt_negative.positions[0]) / 10.0))
        first = ring_a.positions.astype(np.float64)
        second = ring_b.positions.astype(np.float64)
        aromatic_distances.append(float(np.linalg.norm(first.mean(0) - second.mean(0)) / 10.0))
        first_normal = np.cross(first[1] - first[0], first[2] - first[0])
        second_normal = np.cross(second[1] - second[0], second[2] - second[0])
        cosine = abs(np.dot(first_normal, second_normal) / (np.linalg.norm(first_normal) * np.linalg.norm(second_normal)))
        aromatic_angles.append(float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0)))))

    hbond = HydrogenBondAnalysis(
        universe=universe,
        donors_sel="resid 146 and name N",
        hydrogens_sel="resid 146 and (name HN or name HT1 or name HT2 or name HT3)",
        acceptors_sel="resid 151 and name O",
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=150.0,
        update_selections=False,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No hydrogen bonds were found.*")
        hbond.run(verbose=False)
    hbond_frames = {int(frame) for frame in hbond.results.hbonds[:, 0]}

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unlikely unit cell vectors detected.*")
        mdtraj_trajectory = mdtraj.load(str(trajectory_path), top=str(topology_path))
    sasa_atom = mdtraj_trajectory.topology.select("resSeq 197 and name OD2")
    sasa = mdtraj.shrake_rupley(
        mdtraj_trajectory,
        mode="atom",
        probe_radius=0.14,
        n_sphere_points=96,
        change_radii={"H": 0.120, "C": 0.170, "N": 0.155, "O": 0.152, "S": 0.180},
    )[:, sasa_atom[0]]

    checks = [
        event_check(identifier="distance_resid_1_ca_resid_2_ca_lt_0_38_nm", observable="DISTANCE", implementation="NumPy Euclidean norm over MDAnalysis coordinates", query="FIND DISTANCE(resid 1 and name CA, resid 2 and name CA) < 0.38nm;", times=times, states=[value < 0.38 for value in distances], threshold={"comparison": "<", "value": 0.38, "unit": "nm"}),
        event_check(identifier="contact_residue_regions_le_0_45_nm", observable="CONTACT", implementation="MDAnalysis.lib.distances.distance_array minimum", query="FIND CONTACT(resid 30:59, resid 120:159, cutoff=0.45nm);", times=times, states=[value <= 0.45 for value in minimum_contact_distances], threshold={"comparison": "<=", "value": 0.45, "unit": "nm"}),
        event_check(identifier="atom_contact_count_ge_100", observable="CONTACT_COUNT", implementation="NumPy count over an MDAnalysis pair-distance matrix", query="FIND CONTACT_COUNT(resid 30:59, resid 120:159, cutoff=0.45nm, mode=atom) >= 100;", times=times, states=[value >= 100 for value in atom_contact_counts], mode="atom", threshold={"comparison": ">=", "value": 100, "unit": "count"}),
        event_check(identifier="residue_contact_count_ge_5", observable="CONTACT_COUNT", implementation="Python unique residue pairs from an MDAnalysis pair-distance matrix", query="FIND CONTACT_COUNT(resid 30:59, resid 120:159, cutoff=0.45nm, mode=residue) >= 5;", times=times, states=[value >= 5 for value in residue_contact_counts], mode="residue", threshold={"comparison": ">=", "value": 5, "unit": "count"}),
        event_check(identifier="unweighted_rg_protein_lt_1_855_nm", observable="RG", implementation="NumPy geometric-centroid radius of gyration", query="FIND RG(protein) < 1.855nm;", times=times, states=[value < 1.855 for value in unweighted_rg], mass_weighted=False, threshold={"comparison": "<", "value": 1.855, "unit": "nm"}),
        event_check(identifier="mass_weighted_rg_protein_lt_1_855_nm", observable="RG", implementation="NumPy mass-weighted radius of gyration using MDAnalysis masses", query="FIND RG(protein, mass_weighted=true) < 1.855nm;", times=times, states=[value < 1.855 for value in mass_weighted_rg], mass_weighted=True, threshold={"comparison": "<", "value": 1.855, "unit": "nm"}),
        event_check(identifier="hbond_resid_146_n_resid_151_o", observable="HBOND", implementation="MDAnalysis HydrogenBondAnalysis 3.5 A / 150 degree criteria", query="FIND HBOND(resid 146 and name N, resid 151 and name O);", times=times, states=[frame in hbond_frames for frame in range(len(times))], threshold={"distance": 0.35, "distance_unit": "nm", "minimum_angle": 150.0, "angle_unit": "degree"}),
        event_check(identifier="residue_50_phi_lt_minus_90_deg", observable="DIHEDRAL", implementation="MDAnalysis.lib.distances.calc_dihedrals", query="FIND DIHEDRAL(resid 49 and name C, resid 50 and name N, resid 50 and name CA, resid 50 and name C) < -90deg;", times=times, states=[value < -90.0 for value in dihedrals], threshold={"comparison": "<", "value": -90.0, "unit": "degree"}),
        event_check(identifier="protein_helix_fraction_ge_0_54", observable="HELIX", implementation="MDAnalysis dihedrals classified with EnsembleQL's documented phi/psi windows", query="FIND HELIX(protein, minimum_fraction=0.54);", times=times, states=[value >= 0.54 for value in helix_fractions], threshold={"comparison": ">=", "value": 0.54, "unit": "fraction"}),
        event_check(identifier="residue_197_od2_sasa_gt_0_1_nm2", observable="SASA", implementation="MDTraj Shrake-Rupley with 0.14 nm probe and 96 sphere points", query="FIND SASA(resid 197 and name OD2, probe=0.14nm, points=96) > 0.1nm2;", times=times, states=[bool(value > 0.1) for value in sasa], threshold={"comparison": ">", "value": 0.1, "unit": "nm2"}),
        event_check(identifier="fitted_ca_rmsd_lt_0_49_nm", observable="RMSD", implementation="MDAnalysis.analysis.rms.rmsd with centering and superposition", query="FIND RMSD(protein and name CA) < 0.49nm;", times=times, states=[value < 0.49 for value in fitted_rmsds], align=True, threshold={"comparison": "<", "value": 0.49, "unit": "nm"}),
        event_check(identifier="raw_protein_rmsd_lt_0_49_nm", observable="RMSD", implementation="MDAnalysis.analysis.rms.rmsd without centering or superposition", query="FIND RMSD(protein, align=false) < 0.49nm;", times=times, states=[value < 0.49 for value in raw_rmsds], align=False, threshold={"comparison": "<", "value": 0.49, "unit": "nm"}),
        event_check(identifier="coordination_number_ge_0_25", observable="COORDINATION_NUMBER", implementation="NumPy contact count divided by MDAnalysis center selection size", query="FIND COORDINATION_NUMBER(resid 30:59, resid 120:159, cutoff=0.45nm) >= 0.25;", times=times, states=[value >= 0.25 for value in coordination_numbers], threshold={"comparison": ">=", "value": 0.25, "unit": "dimensionless"}),
        event_check(identifier="salt_bridge_resid_78_nh2_resid_79_od1", observable="SALT_BRIDGE", implementation="NumPy Euclidean norm over MDAnalysis coordinates", query="FIND SALT_BRIDGE(resid 78 and name NH2, resid 79 and name OD1);", times=times, states=[value <= 0.4 for value in salt_bridge_distances], threshold={"comparison": "<=", "value": 0.4, "unit": "nm"}),
        event_check(identifier="tyr_181_182_parallel_aromatic_stacking", observable="AROMATIC_STACKING", implementation="NumPy centroid distance and unoriented cross-product normal angle", query="FIND AROMATIC_STACKING(resid 181 and (name CG or name CD1 or name CD2), resid 182 and (name CG or name CD1 or name CD2), distance=0.56nm, max_angle=60deg);", times=times, states=[distance <= 0.56 and angle <= 60.0 for distance, angle in zip(aromatic_distances, aromatic_angles, strict=True)], threshold={"maximum_distance": 0.56, "distance_unit": "nm", "maximum_angle": 60.0, "angle_unit": "degree"}),
    ]

    return {
        "schema_version": 1,
        "dataset": {
            "name": "adenylate kinase DIMS closed-to-open transition",
            "source_package": "MDAnalysisTests",
            "source_version": EXPECTED_VERSION,
            "source_files": {"adk.psf": sha256(Path(PSF)), "adk_dims.dcd": sha256(Path(DCD))},
            "citation_doi": "10.1016/j.jmb.2009.09.009",
            "frame_count": len(times),
            "frame_spacing_ps": 1.0,
            "atom_count": len(universe.atoms),
        },
        "reference_software": {"MDAnalysis": EXPECTED_VERSION, "MDTraj": MDTRAJ_VERSION, "NumPy": np.__version__},
        "prepared_files": {name: sha256(DATA / name) for name in ("topology.pdb", "trajectory.dcd", "trajectory.xtc")},
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if the committed reference is stale")
    args = parser.parse_args()
    generated = build_reference()
    if args.check:
        recorded = json.loads(OUTPUT.read_text())
        if recorded != generated:
            raise SystemExit("The recorded AdK reference is stale; regenerate it without --check")
        print(f"Reference manifest is current: {OUTPUT}")
        return
    OUTPUT.write_text(json.dumps(generated, indent=2) + "\n")
    print(f"Wrote independent reference manifest: {OUTPUT}")


if __name__ == "__main__":
    main()
