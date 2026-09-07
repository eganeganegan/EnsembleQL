from pathlib import Path
from types import SimpleNamespace

import pytest

eql = pytest.importorskip("ensembleql")
pytest.importorskip("ensembleql._core")

ROOT = Path(__file__).parents[2]
EXAMPLE = ROOT / "examples" / "idr_contact_switching"
ADSORPTION = ROOT / "examples" / "peptide_surface_adsorption"
PBC_DATA = ROOT / "tests" / "data"
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


def test_peptide_surface_adsorption_example():
    trajectory = eql.load(
        ADSORPTION / "adsorption.xyz", topology=ADSORPTION / "adsorption.pdb"
    )
    events = trajectory.query((ADSORPTION / "query.eql").read_text())
    assert len(events) == 1
    assert events[0].start == pytest.approx(1.0)
    assert events[0].end == pytest.approx(4.0)

    sasa = trajectory.query("FIND SASA(chain A, points=24) > 0nm2;")
    assert len(sasa) == 1
    assert sasa[0].end == pytest.approx(5.0)


def test_new_temporal_operators_execute():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    immediate = trajectory.query(
        "FIND CONTACT(resid 17, resid 42) "
        "IMMEDIATELY_FOLLOWED_BY CONTACT(resid 17, resid 53);"
    )
    assert len(immediate) == 1
    assert immediate[0].start == pytest.approx(0.0)
    assert immediate[0].end == pytest.approx(4000.0)

    during = trajectory.query(
        "FIND CONTACT(resid 42, resid 53, cutoff=0.01nm) "
        "DURING CONTACT(resid 17, resid 42);"
    )
    assert len(during) == 1
    until = trajectory.query(
        "FIND CONTACT(resid 17, resid 42) UNTIL CONTACT(resid 17, resid 53);"
    )
    assert len(until) == 1
    assert until[0].end == pytest.approx(4000.0)


def test_results_records_and_repr():
    events = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb").query(QUERY)
    assert events.to_records()[0]["duration_ps"] == pytest.approx(4000.0)
    assert "FOLLOWED_BY" in repr(events)


def test_event_analysis_outputs():
    event = lambda start, end, kind: SimpleNamespace(
        start=start, end=end, duration=end - start, type=kind,
        selections=[], metadata={}
    )
    results = eql.EventResults([
        event(0, 1, "A"), event(2, 3, "B"), event(4, 5, "A"),
        event(6, 7, "B"), event(20, 21, "C"),
    ])
    assert results.recurrence_statistics()["A"]["count"] == 2
    assert results.event_frequency(observation_time_ps=10)["A"] == pytest.approx(0.2)
    with pytest.raises(ValueError, match="positive"):
        results.event_frequency(observation_time_ps=0)
    assert results.conditional_probability("A", "B", within_ps=1) == 1.0
    assert results.transition_matrix(within_ps=1)["A"]["B"] == 2
    assert results.motifs(2)[("A", "B")] == 2
    assert results.recurring_subsequences(2) == {("A", "B"): 2}
    assert len(results.temporal_clusters(max_gap_ps=1)) == 2
    assert len(results.event_graph(within_ps=1)["edges"]) == 3
    network = results.state_transition_network(within_ps=1)
    assert {edge["count"] for edge in network["edges"]} == {1, 2}


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


def test_mass_weighted_rg_and_residue_contact_mode():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    assert trajectory.topology.atoms[0].mass_da == pytest.approx(12.011)
    weighted = trajectory.query("FIND RG(protein, mass_weighted=true) < 1nm;")
    assert len(weighted) == 1

    residue_contacts = trajectory.query(
        "FIND CONTACT_COUNT(resid 17, resid 42, mode=residue, cutoff=0.45nm) == 1;"
    )
    assert len(residue_contacts) == 1
    assert residue_contacts[0].end == pytest.approx(2000.0)


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


def test_periodic_minimum_image_query():
    assert eql.minimum_image_distance([0.1, 0.0, 0.0], [1.9, 0.0, 0.0], [2.0, 2.0, 2.0]) == pytest.approx(0.2)
    assert eql.minimum_image_distance_cell(
        [0.1, 0.0, 0.0],
        [1.9, 0.0, 0.0],
        [[2.0, 0.1, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
    ) == pytest.approx(5**0.5 / 10)
    trajectory = eql.load(PBC_DATA / "pbc.xyz", topology=PBC_DATA / "pbc.pdb")
    events = trajectory.query("FIND CONTACT(resid 1, resid 2, cutoff=0.4nm);")
    assert len(events) == 1
    assert events[0].start == pytest.approx(0.0)
    assert events[0].end == pytest.approx(1.0)


def test_periodic_rg_uses_topology_bonds_for_unwrapping():
    trajectory = eql.load(PBC_DATA / "pbc.xyz", topology=PBC_DATA / "pbc.pdb")
    assert trajectory.topology.bonds == [[0, 1]]
    events = trajectory.query("FIND RG(protein) < 0.11nm;")
    assert len(events) == 1
    assert events[0].start == pytest.approx(0.0)
    assert events[0].end == pytest.approx(0.0)


def test_chemfiles_capability_is_discoverable():
    assert isinstance(eql.chemfiles_backend_available(), bool)
    formats = eql.supported_trajectory_extensions()
    assert formats == sorted(formats)
    assert {".ent", ".pdb", ".xyz"} <= set(formats)
    assert (".nc" in formats) is eql.chemfiles_backend_available()


def test_compiled_plan_cache_is_trajectory_scoped():
    trajectory = eql.load(EXAMPLE / "switching.xyz", topology=EXAMPLE / "switching.pdb")
    query = "FIND CONTACT(resid 17, resid 42);"
    assert trajectory.cached_plan_count == 0
    trajectory.query(query)
    trajectory.query(query)
    trajectory.explain(query)
    assert trajectory.cached_plan_count == 1
    trajectory.query("FIND CONTACT(resid 17, resid 53);")
    assert trajectory.cached_plan_count == 2
    trajectory.clear_plan_cache()
    assert trajectory.cached_plan_count == 0


def test_configurable_default_timestep():
    trajectory = eql.load(
        PBC_DATA / "no_time.xyz",
        topology=PBC_DATA / "pbc.pdb",
        default_timestep="2.5ps",
    )
    events = trajectory.query("FIND CONTACT(resid 1, resid 2, cutoff=0.4nm);")
    assert len(events) == 1
    assert events[0].end == pytest.approx(2.5)

    with pytest.raises((ValueError, RuntimeError), match="positive"):
        eql.load(
            PBC_DATA / "no_time.xyz",
            topology=PBC_DATA / "pbc.pdb",
            default_timestep="0ps",
        )


def test_multimodel_pdb_trajectory():
    trajectory = eql.load(
        PBC_DATA / "models.pdb",
        topology=PBC_DATA / "models.pdb",
        default_timestep="2.5ps",
    )
    assert len(trajectory.topology.atoms) == 2
    assert trajectory.topology.bonds == [[0, 1]]
    events = trajectory.query("FIND CONTACT(resid 1, resid 2, cutoff=0.4nm);")
    assert len(events) == 1
    assert events[0].start == pytest.approx(0.0)
    assert events[0].end == pytest.approx(2.5)
