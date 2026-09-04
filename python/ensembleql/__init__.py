"""Python interface for the EnsembleQL temporal molecular query engine."""

from .trajectory import Trajectory, load
from .results import EventResults

try:
    from ._core import Atom, Event, QueryError, Topology
except ImportError:  # A source checkout may not have built the native module yet.
    Atom = Event = Topology = None

    class QueryError(RuntimeError):
        pass

__all__ = ["Atom", "Event", "EventResults", "QueryError", "Topology", "Trajectory", "load"]
__version__ = "0.1.0"
