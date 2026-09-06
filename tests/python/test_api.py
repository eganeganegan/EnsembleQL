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


def test_boolean_selection_and_custom_cutoff():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    assert trajectory.topology.select("protein and name CA") == [0, 1, 2]
    events = trajectory.query("FIND CONTACT(resid 17, resid 42, cutoff=0.35nm);")
    assert len(events) == 1
    assert events[0].end == pytest.approx(1000.0)


def test_cutoff_rejects_time_dimension():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    with pytest.raises(eql.QueryError, match="Distance expected, received time unit"):
        trajectory.query("FIND CONTACT(resid 17, resid 42, cutoff=2ns);")


def test_explain_returns_structured_streaming_plan():
    plan = eql.explain(QUERY, topology=EXAMPLE / "switching.pdb")
    assert plan["streaming"] is True
    assert len(plan["selections"]) == 3
    assert len(plan["observables"]) == 2
    assert plan["temporal_operations"] == ["FOLLOWED_BY WITHIN 3000ps"]
    assert "FOLLOWED_BY" in plan["plan_tree"]

    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    assert trajectory.explain(QUERY) == plan
