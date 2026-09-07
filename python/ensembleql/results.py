"""Result containers independent of optional dataframe libraries."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from math import isfinite
from statistics import fmean, median


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

    def recurrence_statistics(self) -> dict[str, dict[str, float | int]]:
        """Summarize event counts, durations, and inter-event gaps by type."""
        grouped: dict[str, list] = {}
        for event in self._events:
            grouped.setdefault(event.type, []).append(event)
        result = {}
        for event_type, events in grouped.items():
            ordered = sorted(events, key=lambda event: (event.start, event.end))
            durations = [event.duration for event in ordered]
            gaps = [
                max(0.0, current.start - previous.end)
                for previous, current in zip(ordered, ordered[1:])
            ]
            result[event_type] = {
                "count": len(ordered),
                "total_duration_ps": sum(durations),
                "mean_duration_ps": fmean(durations),
                "median_duration_ps": median(durations),
                "mean_gap_ps": fmean(gaps) if gaps else 0.0,
            }
        return result

    def event_frequency(self, *, observation_time_ps: float) -> dict[str, float]:
        """Return occurrence rates by event type in events per picosecond."""
        duration = float(observation_time_ps)
        if not isfinite(duration) or duration <= 0:
            raise ValueError("observation time must be finite and positive")
        counts: dict[str, int] = {}
        for event in self._events:
            counts[event.type] = counts.get(event.type, 0) + 1
        return {event_type: count / duration for event_type, count in counts.items()}

    def conditional_probability(
        self, antecedent: str, consequent: str, *, within_ps: float
    ) -> float:
        """Estimate P(consequent follows antecedent within an inclusive window)."""
        window = self._validated_window(within_ps)
        left = [event for event in self._events if event.type == antecedent]
        right = [event for event in self._events if event.type == consequent]
        if not left:
            return 0.0
        matches = sum(
            any(
                event.end <= candidate.start <= event.end + window
                for candidate in right
            )
            for event in left
        )
        return matches / len(left)

    def transition_matrix(self, *, within_ps: float | None = None) -> dict[str, dict[str, int]]:
        """Count transitions between chronologically adjacent event types."""
        window = None if within_ps is None else self._validated_window(within_ps)
        ordered = self._ordered_events()
        event_types = sorted({event.type for event in ordered})
        matrix = {source: {target: 0 for target in event_types} for source in event_types}
        for source, target in zip(ordered, ordered[1:]):
            gap = max(0.0, target.start - source.end)
            if window is None or gap <= window:
                matrix[source.type][target.type] += 1
        return matrix

    def motifs(
        self, length: int, *, within_ps: float | None = None
    ) -> dict[tuple[str, ...], int]:
        """Count sliding event-type subsequences of a fixed positive length."""
        if not isinstance(length, int) or isinstance(length, bool) or length <= 0:
            raise ValueError("motif length must be a positive integer")
        window = None if within_ps is None else self._validated_window(within_ps)
        ordered = self._ordered_events()
        result: dict[tuple[str, ...], int] = {}
        for first in range(len(ordered) - length + 1):
            group = ordered[first : first + length]
            if window is not None and group[-1].end - group[0].start > window:
                continue
            motif = tuple(event.type for event in group)
            result[motif] = result.get(motif, 0) + 1
        return result

    def recurring_subsequences(
        self, length: int, *, minimum_count: int = 2, within_ps: float | None = None
    ) -> dict[tuple[str, ...], int]:
        """Return motifs occurring at least ``minimum_count`` times."""
        if not isinstance(minimum_count, int) or isinstance(minimum_count, bool) or minimum_count <= 0:
            raise ValueError("minimum_count must be a positive integer")
        return {
            motif: count
            for motif, count in self.motifs(length, within_ps=within_ps).items()
            if count >= minimum_count
        }

    def temporal_clusters(self, *, max_gap_ps: float) -> list[dict]:
        """Group events connected by gaps no larger than ``max_gap_ps``."""
        maximum_gap = self._validated_window(max_gap_ps)
        ordered = self._ordered_events()
        if not ordered:
            return []
        clusters = []
        members = [0]
        cluster_end = ordered[0].end
        for index, event in enumerate(ordered[1:], start=1):
            if event.start - cluster_end <= maximum_gap:
                members.append(index)
                cluster_end = max(cluster_end, event.end)
                continue
            clusters.append(self._cluster_record(ordered, members))
            members = [index]
            cluster_end = event.end
        clusters.append(self._cluster_record(ordered, members))
        return clusters

    def event_graph(self, *, within_ps: float | None = None) -> dict[str, list[dict]]:
        """Create a directed graph of temporally ordered events."""
        window = None if within_ps is None else self._validated_window(within_ps)
        ordered = self._ordered_events()
        nodes = [
            {"id": index, "type": event.type, "start_ps": event.start, "end_ps": event.end}
            for index, event in enumerate(ordered)
        ]
        edges = []
        for source, event in enumerate(ordered):
            for target in range(source + 1, len(ordered)):
                candidate = ordered[target]
                if candidate.start < event.end:
                    continue
                gap = candidate.start - event.end
                if window is not None and gap > window:
                    break
                edges.append({"source": source, "target": target, "gap_ps": gap})
        return {"nodes": nodes, "edges": edges}

    def state_transition_network(self, *, within_ps: float | None = None) -> dict[str, list[dict]]:
        """Aggregate adjacent event transitions into a weighted state graph."""
        matrix = self.transition_matrix(within_ps=within_ps)
        nodes = [{"id": event_type} for event_type in matrix]
        edges = [
            {"source": source, "target": target, "count": count}
            for source, targets in matrix.items()
            for target, count in targets.items()
            if count
        ]
        return {"nodes": nodes, "edges": edges}

    def _ordered_events(self) -> list:
        return sorted(self._events, key=lambda event: (event.start, event.end, event.type))

    @staticmethod
    def _validated_window(value: float) -> float:
        value = float(value)
        if not isfinite(value) or value < 0:
            raise ValueError("time window must be finite and non-negative")
        return value

    @staticmethod
    def _cluster_record(events: list, members: list[int]) -> dict:
        selected = [events[index] for index in members]
        return {
            "start_ps": min(event.start for event in selected),
            "end_ps": max(event.end for event in selected),
            "event_indices": members,
            "types": [event.type for event in selected],
        }

    def __repr__(self) -> str:
        if not self._events:
            return "EnsembleQL events: 0"
        lines = ["type          start (ps)   end (ps)   duration (ps)"]
        for event in self._events:
            lines.append(f"{event.type:<13} {event.start:>10.3f} {event.end:>10.3f} {event.duration:>15.3f}")
        return "\n".join(lines)
