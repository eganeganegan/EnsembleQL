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
