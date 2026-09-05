# Architecture notes

`FrameReader` owns incremental decoding; `Trajectory` combines it with immutable topology metadata. Execution resets the reader, validates strictly increasing timestamps, evaluates predicate sources once per frame, and updates an `EventExtractor` state machine per source. A per-frame value cache is keyed by canonical observable identity, so repeated and symmetrically reversed observables are computed once even when multiple comparisons consume them. Only completed interval events are retained.

The parser emits typed AST nodes and performs dimensional validation. The planner resolves every required selection early, reports zero-match errors before trajectory scanning, and records unique observables. Temporal evaluation is a separate algebra over `Event` values, keeping molecular geometry unaware of relations such as `FOLLOWED_BY`.

Future geometry kernels should replace implementations, not interfaces. Future compressed formats should implement `FrameReader`; they must expose box information so unsupported PBC can never silently degrade into Euclidean non-PBC analysis.
