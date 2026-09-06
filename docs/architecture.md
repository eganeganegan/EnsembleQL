# Architecture notes

`FrameReader` owns incremental decoding; `Trajectory` combines it with immutable topology metadata. Execution resets the reader, validates strictly increasing timestamps, evaluates predicate sources once per frame, and updates an `EventExtractor` state machine per source. A per-frame value cache is keyed by canonical observable identity, so repeated and symmetrically reversed observables are computed once even when multiple comparisons consume them. Only completed interval events are retained.

The parser emits typed AST nodes and performs dimensional validation. The planner resolves every required selection early, reports zero-match errors before trajectory scanning, and records unique observables. Temporal evaluation is a separate algebra over `Event` values, keeping molecular geometry unaware of relations such as `FOLLOWED_BY`.

Future geometry kernels should replace implementations, not interfaces. Orthorhombic frame boxes currently select minimum-image distance/contact geometry automatically. Triclinic boxes and observables requiring molecule unwrapping are rejected explicitly. Future compressed formats should implement `FrameReader` and must always expose box information so periodic data can never silently degrade into non-periodic geometry.
