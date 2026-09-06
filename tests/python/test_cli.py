import json
from pathlib import Path

import pytest

pytest.importorskip("ensembleql._core")
from ensembleql.cli import main

ROOT = Path(__file__).parents[2]
EXAMPLE = ROOT / "examples" / "idr_contact_switching"


def test_cli_json(capsys):
    status = main([
        "query", "--topology", str(EXAMPLE / "switching.pdb"),
        "--trajectory", str(EXAMPLE / "switching.xyz"),
        "--file", str(EXAMPLE / "query.eql"), "--format", "json",
    ])
    assert status == 0
    assert json.loads(capsys.readouterr().out)[0]["type"] == "FOLLOWED_BY"


def test_cli_explain_json(capsys):
    status = main([
        "explain", "--topology", str(EXAMPLE / "switching.pdb"),
        "--file", str(EXAMPLE / "query.eql"), "--format", "json",
    ])
    assert status == 0
    plan = json.loads(capsys.readouterr().out)
    assert plan["streaming"] is True
    assert plan["temporal_operations"] == ["FOLLOWED_BY WITHIN 3000ps"]


def test_cli_explain_text(capsys):
    status = main([
        "explain", "--topology", str(EXAMPLE / "switching.pdb"),
        "--query", "FIND CONTACT(resid 17, resid 42);",
    ])
    assert status == 0
    output = capsys.readouterr().out
    assert "Streaming: yes" in output
    assert "Plan:" in output
    assert "CONTACT" in output


def test_cli_default_timestep(capsys):
    data = ROOT / "tests" / "data"
    status = main([
        "query", "--topology", str(data / "pbc.pdb"),
        "--trajectory", str(data / "no_time.xyz"),
        "--default-timestep", "2.5ps",
        "--query", "FIND CONTACT(resid 1, resid 2, cutoff=0.4nm);",
        "--format", "json",
    ])
    assert status == 0
    assert json.loads(capsys.readouterr().out)[0]["end_ps"] == pytest.approx(2.5)
