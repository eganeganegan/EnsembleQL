import hashlib
import json
from pathlib import Path

import pytest


eql = pytest.importorskip("ensembleql")
pytest.importorskip("ensembleql._core")

ROOT = Path(__file__).parents[2]
FIXTURES = {
    "adk": {
        "data": ROOT / "validation" / "data" / "adk",
        "reference": json.loads(
            (ROOT / "validation" / "adk_reference.json").read_text()
        ),
        "prepare": "validation/prepare_adk.py",
    },
    "membrane-peptide": {
        "data": ROOT / "validation" / "data" / "membrane-peptide",
        "reference": json.loads(
            (ROOT / "validation" / "membrane_peptide_reference.json").read_text()
        ),
        "prepare": "validation/prepare_membrane_peptide.py --download",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_fixture(name: str) -> tuple[Path, dict]:
    if not eql.chemfiles_backend_available():
        pytest.skip("scientific validation requires the Chemfiles backend")

    fixture = FIXTURES[name]
    data = fixture["data"]
    reference = fixture["reference"]
    missing = [
        filename
        for filename in reference["prepared_files"]
        if not (data / filename).is_file()
    ]
    if missing:
        pytest.skip(
            f"prepare the external {name} fixture with {fixture['prepare']}; "
            f"missing: {', '.join(missing)}"
        )

    for filename, expected in reference["prepared_files"].items():
        assert sha256(data / filename) == expected, (
            f"unexpected validation fixture: {filename}"
        )
    return data, reference


CASES = [
    (fixture_name, check)
    for fixture_name, fixture in FIXTURES.items()
    for check in fixture["reference"]["checks"]
]


@pytest.mark.parametrize(
    ("fixture_name", "check"),
    CASES,
    ids=[f"{fixture}-{check['id']}" for fixture, check in CASES],
)
def test_events_match_independent_reference(fixture_name, check):
    data, reference = require_fixture(fixture_name)
    trajectory = eql.load(data / "trajectory.xtc", topology=data / "topology.pdb")

    assert len(trajectory.topology.atoms) == reference["dataset"]["atom_count"]
    events = trajectory.query(check["query"])
    assert len(events) == check["event_count"]

    actual = [(event.start, event.end) for event in events]
    expected = [tuple(interval) for interval in check["intervals_ps"]]
    tolerance = check["time_tolerance_ps"]
    for actual_interval, expected_interval in zip(actual, expected, strict=True):
        assert actual_interval[0] == pytest.approx(expected_interval[0], abs=tolerance)
        assert actual_interval[1] == pytest.approx(expected_interval[1], abs=tolerance)

    frame_spacing = reference["dataset"]["frame_spacing_ps"]
    matching_frames = sum(
        round((end - start) / frame_spacing) + 1 for start, end in actual
    )
    assert matching_frames == check["matching_frame_count"]


def test_reference_manifests_cover_every_observable_family():
    covered = {
        check["observable"]
        for fixture in FIXTURES.values()
        for check in fixture["reference"]["checks"]
    }
    assert covered == {
        "AROMATIC_STACKING",
        "CONTACT",
        "CONTACT_COUNT",
        "COORDINATION_NUMBER",
        "DIHEDRAL",
        "DISTANCE",
        "HBOND",
        "HELIX",
        "ORIENTATION",
        "RG",
        "RMSD",
        "SALT_BRIDGE",
        "SASA",
        "SURFACE_DISTANCE",
    }
