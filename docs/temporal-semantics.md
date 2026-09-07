# Sampled-time event semantics

EnsembleQL reports what is observed at trajectory samples. It does not interpolate molecular states between frames.

## Frame times

Frame times are finite values in ps and must be strictly increasing. Irregular spacing is supported. Duplicate or decreasing timestamps are rejected before a result is returned. XYZ comments may include an explicit value such as `time=2.5ns`; without one, the reader assigns the configured fallback step from frame zero. The fallback defaults to 1 ps and must be finite and positive. Python and CLI callers provide it as an explicit-unit duration such as `2fs` or `0.5ns`. File-provided timestamps take precedence frame by frame.

## Event intervals

A frame predicate produces an event for each maximal contiguous run of true samples. An event starts at the first true sample and ends at the last true sample. Its observed duration is:

```text
duration = last_true_sample_time - first_true_sample_time
```

Consequently, a one-sample event has zero observed duration. EnsembleQL does not extend it to the next frame or infer how long it remained true after sampling. On irregularly sampled data, `FOR >= 10ps` compares 10 ps against this observed span—not against a frame count.

## Boundaries

All current duration, cutoff, and temporal window boundaries are inclusive:

- a distance exactly at the contact cutoff is a contact;
- an event lasting exactly the requested `FOR` duration passes;
- a transition gap exactly equal to `WITHIN` passes;
- closed intervals sharing an endpoint both overlap and satisfy non-strict `BEFORE`/`AFTER` ordering.

The geometry implementation allows only a tiny floating-point tolerance at inclusive distance boundaries to prevent angstrom-to-nm conversion roundoff from changing classification.

`FOLLOWED_BY` pairs each left event with the earliest right event whose start is at or after the left event's end and whose non-negative gap satisfies `WITHIN`, when present. A shared sampled endpoint therefore has a zero-ps gap.

Unparenthesized temporal operators chain left-associatively. Thus `A FOLLOWED_BY B WITHIN 1ps FOLLOWED_BY C WITHIN 2ps` first constructs the `A`→`B` events and then joins those to `C`. Parentheses can make any alternative grouping explicit. Frame-level `AND` binds more tightly than temporal operators.

Additional interval operators have the following definitions:

- `DURING` returns pairs where the closed left interval is contained in the closed right interval, including equal endpoints.
- `PRECEDES` is the strict counterpart to `BEFORE`: the left event must end before the right event starts, so touching intervals do not qualify.
- `IMMEDIATELY_FOLLOWED_BY` requires the left end and right start to be the same sampled timestamp, within floating-point time tolerance.
- `UNTIL` pairs each left event with the earliest right event beginning at or after the left event ends. It describes observed event ordering and does not infer an unobserved state in any sampling gap.
- `REPEATS >= N` emits each sliding window of `N` consecutive source events. `WITHIN` optionally limits the inclusive span from the first event's start through the last event's end.

These definitions are deliberately explicit and covered by native golden tests. Alternative interval estimators may be added later, but must be opt-in and labeled in result metadata.
