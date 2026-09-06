"""Trajectory loading and query convenience API."""

from __future__ import annotations

from pathlib import Path

from .results import EventResults


class Trajectory:
    """A streaming native trajectory and its molecular topology."""

    def __init__(self, trajectory: str | Path, topology: str | Path):
        try:
            from ._core import NativeTrajectory
        except ImportError as exc:
            raise RuntimeError(
                "EnsembleQL's native extension is not built. Install the project with "
                "`python -m pip install -e .` before loading trajectories."
            ) from exc
        self._native = NativeTrajectory.from_files(str(trajectory), str(topology))

    @property
    def topology(self):
        return self._native.topology

    def query(self, text: str) -> EventResults:
        """Parse, plan, and execute *text* in one streaming pass."""
        return EventResults(self._native.query(text))

    def explain(self, text: str) -> dict:
        """Return the resolved execution plan without reading trajectory frames."""
        return dict(self._native.explain(text))


def load(trajectory: str | Path, *, topology: str | Path) -> Trajectory:
    """Open an XYZ trajectory using atom metadata from a PDB topology."""
    return Trajectory(trajectory, topology)
