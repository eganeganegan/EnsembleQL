"""Python interface for the EnsembleQL temporal molecular query engine."""

from .trajectory import Trajectory, load
from .results import EventResults
from .query import explain

try:
    from ._core import Atom, Event, QueryError, Topology, minimum_image_distance
except ImportError:  # A source checkout may not have built the native module yet.
    Atom = Event = Topology = minimum_image_distance = None

    class QueryError(RuntimeError):
        pass

__all__ = ["Atom", "Event", "EventResults", "QueryError", "Topology", "Trajectory", "explain", "load", "minimum_image_distance"]
__version__ = "0.1.0"
