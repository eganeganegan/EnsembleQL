# Architecture notes

`FrameReader` owns incremental decoding; `Trajectory` combines it with immutable topology metadata. Execution resets the reader, evaluates deduplicated predicate sources once per frame, and updates an `EventExtractor` state machine per source. Only completed interval events are retained.

The parser emits typed AST nodes and performs dimensional validation. The planner resolves every required selection early, reports zero-match errors before trajectory scanning, and records unique observables. Temporal evaluation is a separate algebra over `Event` values, keeping molecular geometry unaware of relations such as `FOLLOWED_BY`.

Future geometry kernels should replace implementations, not interfaces. Future compressed formats should implement `FrameReader`; they must expose box information so unsupported PBC can never silently degrade into Euclidean non-PBC analysis.
