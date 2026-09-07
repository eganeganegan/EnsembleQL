"""Trajectory loading and query convenience API."""

from __future__ import annotations

from pathlib import Path

from .results import EventResults


class Trajectory:
    """A streaming native trajectory and its molecular topology."""

    def __init__(
        self,
        trajectory: str | Path,
        topology: str | Path,
        *,
        default_timestep: str = "1ps",
    ):
        try:
            from ._core import NativeTrajectory, parse_duration_ps
        except ImportError as exc:
            raise RuntimeError(
                "EnsembleQL's native extension is not built. Install the project with "
                "`python -m pip install -e .` before loading trajectories."
            ) from exc
        if not isinstance(default_timestep, str):
            raise TypeError("default_timestep must be a duration string with explicit units")
        timestep_ps = parse_duration_ps(default_timestep)
        if timestep_ps <= 0:
            raise ValueError("default_timestep must be positive")
        self._native = NativeTrajectory.from_files(
            str(trajectory), str(topology), timestep_ps
        )

    @property
    def topology(self):
        return self._native.topology

    def query(self, text: str) -> EventResults:
        """Compile or reuse *text*, then execute it in one streaming pass."""
        return EventResults(self._native.query(text))

    def explain(self, text: str) -> dict:
        """Compile or reuse *text* and return its plan without reading frames."""
        return dict(self._native.explain(text))

    @property
    def cached_plan_count(self) -> int:
        """Number of compiled query plans retained for this topology."""
        return self._native.cached_plan_count

    def clear_plan_cache(self) -> None:
        """Discard all compiled query plans retained by this trajectory."""
        self._native.clear_plan_cache()


def load(
    trajectory: str | Path,
    *,
    topology: str | Path,
    default_timestep: str = "1ps",
) -> Trajectory:
    """Open a supported trajectory using atom metadata and bonds from a PDB topology.

    Use ``ensembleql.supported_trajectory_extensions()`` to inspect the formats
    enabled in the installed native extension.
    """
    return Trajectory(trajectory, topology, default_timestep=default_timestep)
