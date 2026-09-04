# EnsembleQL

EnsembleQL is an open-source, declarative temporal query engine for molecular-dynamics trajectories. It lets computational biologists and chemists ask when molecular behavior occurs and how events relate in time, instead of rebuilding each analysis as a bespoke array-processing script.

```text
trajectory -> observables -> predicates -> events -> temporal relationships
```

The first milestone focuses on transient contacts in intrinsically disordered proteins (IDRs). Its execution engine is C++20; pybind11 exposes a compact Python and CLI interface.

## Why events, not only averages?

Frame-wise contact probabilities can erase order. A conventional analysis might say:

```text
R17-D42 = 32%
R17-E53 = 29%
```

Those values cannot distinguish independent contacts from a directed interaction switch. EnsembleQL can instead report:

```text
R17-D42 -> R17-E53 switching events: 14
median transition gap: 0.8 ns
```

Events are maximal contiguous intervals over which a frame predicate is true. Temporal operators act on those intervals independently of the molecular observable that produced them.

## Quick start

Build and test the native core:

```bash
cmake -S . -B build -DENSEMBLEQL_BUILD_PYTHON=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For the Python package (Python 3.11+), install into an isolated environment. The build installs its pybind11/scikit-build dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
pytest
```

Run the included IDR contact-switching example:

```python
from pathlib import Path
import ensembleql as eql

root = Path("examples/idr_contact_switching")
traj = eql.load(root / "switching.xyz", topology=root / "switching.pdb")
events = traj.query("""
    FIND CONTACT(resid 17, resid 42)
    FOLLOWED_BY CONTACT(resid 17, resid 53)
    WITHIN 3ns;
""")
print(events)
```

The synthetic trajectory includes an R17-D42 contact, a one-frame boundary where both contacts exist, then an R17-E53 contact. EnsembleQL returns one directed switch from 0 to 4000 ps with a zero-ps transition gap. This preserves timing and direction that two contact probabilities do not.

## DSL

Implemented queries include:

```text
FIND CONTACT(resid 17, resid 42);
FIND CONTACT(resid 17, resid 42) FOR >= 5ns;
FIND DISTANCE(resid 17, resid 42) < 0.8nm;
FIND RG(protein) < 2.0nm FOR >= 10ns;
FIND CONTACT_COUNT(resid 1:20, resid 40:60) >= 4;
FIND CONTACT(resid 17, resid 42) OVERLAPS CONTACT(resid 17, resid 53);
```

`FOLLOWED_BY`, `WITHIN`, `BEFORE`, `AFTER`, `OVERLAPS`, `FOR`, `AND`, and `OR` have AST nodes. Parentheses can group temporal expressions. The parser only creates an AST; the planner resolves selections and deduplicates required observables; the engine then evaluates every required predicate in a single streaming traversal and retains events rather than a full boolean time series.

Selections currently support `resid 17`, `resid 17:25`, `name CA`, `resname ARG`, `chain A`, and `protein`. The internal interfaces leave boolean selection parsing as a later extension.

## CLI

```bash
ensembleql query \
  --topology examples/idr_contact_switching/switching.pdb \
  --trajectory examples/idr_contact_switching/switching.xyz \
  --file examples/idr_contact_switching/query.eql

ensembleql query --topology structure.pdb --trajectory trajectory.xyz \
  --query "FIND CONTACT(resid 17, resid 42) FOR >= 2ns;" --format json
```

Output formats are `table`, `csv`, and `json`. `EventResults.to_dataframe()` returns a pandas DataFrame when pandas is installed and otherwise returns a list of records.

## Scientific definitions and units

- Coordinates and distances are normalized to nm. XYZ coordinates are interpreted as angstroms; PDB is used for topology metadata only.
- Times are normalized to ps. Supported distance units are `nm`, `angstrom`, and `A`; supported time units are `fs`, `ps`, `ns`, and `us`. Query thresholds require explicit units except integer-like counts.
- `DISTANCE(A,B)` is the minimum Euclidean distance between distinct atoms in A and B.
- `CONTACT(A,B)` is true when that minimum distance is less than or equal to 0.45 nm. `CONTACT_COUNT` counts atom pairs at or below the same inclusive cutoff.
- `RG(A)` is the unweighted root-mean-square distance of selected atom coordinates from their geometric centroid. Mass-weighting is not yet implemented.
- Event endpoints are the timestamps of the first and last true sampled frames; duration is `end - start`. `FOR >=`, contact cutoffs, `WITHIN`, and interval boundary comparisons are inclusive.
- Periodic boundary conditions are unsupported. A frame declaring a box is rejected rather than silently analyzed with non-periodic distances.

## Architecture

Public headers separate trajectory/topology I/O, selections, geometry, observables, event extraction, interval algebra, AST parsing, planning, and execution. `FrameReader` is the backend-neutral streaming interface; future XTC/TRR/DCD readers can implement it without changing query execution. Observable classes are likewise independent of the parser.

The straightforward contact kernel is currently O(N×M). Its stable API permits cell lists, spatial hashing, neighbor lists, SIMD, or OpenMP underneath it. Other intended extension points are parallel frame evaluation, thread pools, memory-mapped or compressed trajectory readers, and query-plan caching. See [ROADMAP.md](ROADMAP.md).

## Current limitations

Only PDB topology metadata and XYZ trajectories are supported. There is no PBC, mass-weighted RG, configurable contact cutoff syntax, full boolean selection grammar, or compressed trajectory reader yet. Event interval semantics use sampled timestamps and therefore do not infer behavior between frames. `AND`/`OR` combine frame predicates; temporal relations operate on extracted intervals.

## Development

The native target compiles with `-Wall -Wextra -Wpedantic`. Enable microbenchmarks with `-DENSEMBLEQL_BUILD_BENCHMARKS=ON`. The synthetic benchmarks cover contact detection at increasing atom counts, streaming event extraction, and temporal joins.

Contributions should preserve scientific definitions, add boundary-condition tests, and keep file-format backends independent from the engine.
