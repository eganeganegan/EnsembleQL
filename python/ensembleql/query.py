"""Query validation helpers."""

from pathlib import Path


def validate(text: str) -> None:
    """Raise QueryError if *text* is not valid EnsembleQL syntax."""
    try:
        from ._core import parse_query
    except ImportError as exc:
        raise RuntimeError("EnsembleQL's native extension is not built") from exc
    parse_query(text)


def explain(text: str, *, topology: str | Path) -> dict:
    """Plan *text* against a PDB topology without opening a trajectory."""
    try:
        from ._core import Topology, explain_query
    except ImportError as exc:
        raise RuntimeError("EnsembleQL's native extension is not built") from exc
    return dict(explain_query(Topology.from_pdb(str(topology)), text))
