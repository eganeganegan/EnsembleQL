"""Result containers independent of optional dataframe libraries."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence


class EventResults(Sequence):
    def __init__(self, events: Iterable):
        self._events = list(events)

    def __len__(self) -> int:
        return len(self._events)

    def __getitem__(self, index):
        return self._events[index]

    def __iter__(self) -> Iterator:
        return iter(self._events)

    def to_records(self) -> list[dict]:
        return [
            {
                "start_ps": event.start,
                "end_ps": event.end,
                "duration_ps": event.duration,
                "type": event.type,
                "selections": list(event.selections),
                "metadata": dict(event.metadata),
            }
            for event in self._events
        ]

    def to_dataframe(self):
        """Return pandas.DataFrame when pandas is installed, otherwise records."""
        records = self.to_records()
        try:
            import pandas as pd
        except ImportError:
            return records
        return pd.DataFrame.from_records(records)

    def __repr__(self) -> str:
        if not self._events:
            return "EnsembleQL events: 0"
        lines = ["type          start (ps)   end (ps)   duration (ps)"]
        for event in self._events:
            lines.append(f"{event.type:<13} {event.start:>10.3f} {event.end:>10.3f} {event.duration:>15.3f}")
        return "\n".join(lines)
