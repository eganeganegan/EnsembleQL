from pathlib import Path

import pytest

eql = pytest.importorskip("ensembleql")
pytest.importorskip("ensembleql._core")

ROOT = Path(__file__).parents[2]
EXAMPLE = ROOT / "examples" / "idr_contact_switching"
QUERY = """
FIND CONTACT(resid 17, resid 42)
FOLLOWED_BY CONTACT(resid 17, resid 53)
WITHIN 3ns;
"""


def test_contact_switching_vertical_slice():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    events = trajectory.query(QUERY)
    assert len(events) == 1
    assert events[0].start == pytest.approx(0.0)
    assert events[0].end == pytest.approx(4000.0)
    assert events[0].metadata["transition_gap_ps"] == "0.000000"


def test_results_records_and_repr():
    events = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb").query(QUERY)
    assert events.to_records()[0]["duration_ps"] == pytest.approx(4000.0)
    assert "FOLLOWED_BY" in repr(events)


def test_engine_error_reaches_python():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    with pytest.raises(RuntimeError, match="matched zero atoms"):
        trajectory.query("FIND CONTACT(resid 17, resid 999);")
