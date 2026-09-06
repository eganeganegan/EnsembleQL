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

Build and test the dependency-free native core:

```bash
cmake -S . -B build -DENSEMBLEQL_BUILD_PYTHON=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Enable XTC, TRR, and DCD through an installed chemfiles library, or fetch the pinned stable release during configuration:

```bash
cmake -S . -B build -DENSEMBLEQL_FETCH_CHEMFILES=ON
cmake --build build -j
```

Use `-DENSEMBLEQL_REQUIRE_CHEMFILES=ON` when configuration should fail rather than produce an XYZ-only build. Python exposes `eql.chemfiles_backend_available()` for capability checks.

To include the backend in an editable Python installation:

```bash
CMAKE_ARGS="-DENSEMBLEQL_FETCH_CHEMFILES=ON -DENSEMBLEQL_REQUIRE_CHEMFILES=ON" \
  python -m pip install -e .
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
FIND CONTACT(resid 17, resid 42, cutoff=0.35nm);
FIND CONTACT(resid 17, resid 42) FOR >= 5ns;
FIND DISTANCE(resid 17, resid 42) < 0.8nm;
FIND RG(protein) < 2.0nm FOR >= 10ns;
FIND CONTACT_COUNT(resid 1:20, resid 40:60) >= 4;
FIND CONTACT(resid 17, resid 42) OVERLAPS CONTACT(resid 17, resid 53);
```

`FOLLOWED_BY`, `WITHIN`, `BEFORE`, `AFTER`, `OVERLAPS`, `FOR`, `AND`, and `OR` have AST nodes. Parentheses can group temporal expressions. The parser only creates an AST; the planner resolves selections and deduplicates required observables; the engine then evaluates every required predicate in a single streaming traversal and retains events rather than a full boolean time series.

Selections support `resid 17`, `resid 17:25`, `name CA`, `resname ARG`, `chain A`, `protein`, and `hydrogen`. They compose with case-insensitive `and`, `or`, `not`, and parentheses; `and` binds more tightly than `or`.

`CONTACT` and `CONTACT_COUNT` use a 0.45 nm default cutoff. Override it with a dimensionally checked third argument such as `cutoff=4A` or `cutoff=0.35nm`.

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

Inspect a query before reading trajectory frames:

```bash
ensembleql explain \
  --topology examples/idr_contact_switching/switching.pdb \
  --file examples/idr_contact_switching/query.eql
```

The explanation reports resolved selection expressions, canonical deduplicated observables, frame predicates, temporal operations, and the execution-plan tree. Use `--format json` for machine-readable output. Python provides the same information through `traj.explain(query)` or `eql.explain(query, topology="structure.pdb")`.

## Scientific definitions and units

- Coordinates and distances are normalized to nm. XYZ coordinates are interpreted as angstroms; PDB is used for topology metadata only.
- Times are normalized to ps. Supported distance units are `nm`, `angstrom`, and `A`; supported time units are `fs`, `ps`, `ns`, and `us`. Query thresholds require explicit units except integer-like counts.
- `DISTANCE(A,B)` is the minimum distance between distinct atoms in A and B. It uses Euclidean distance without a cell and a nearest-image search for orthorhombic or triclinic periodic cells.
- `CONTACT(A,B)` is true when that minimum distance is less than or equal to 0.45 nm. `CONTACT_COUNT` counts atom pairs at or below the same inclusive cutoff. Both honor orthorhombic and triclinic periodic boundaries.
- `RG(A)` is the unweighted root-mean-square distance of selected atom coordinates from their geometric centroid. On periodic frames, the selected atoms are reconstructed by traversing PDB `CONECT` bonds before calculating the centroid. Mass-weighting is not yet implemented.
- Event endpoints are the timestamps of the first and last true sampled frames; duration is `end - start`. Single-sample events therefore have zero observed duration. Irregular timestamps are supported, while duplicate or decreasing timestamps are rejected. `FOR >=`, contact cutoffs, `WITHIN`, and interval boundary comparisons are inclusive. See [sampled-time event semantics](docs/temporal-semantics.md).
- Periodic cells use a nearest-image search over neighboring lattice translations. XYZ comments accept `box=20,20,20A`; extended XYZ accepts a full `Lattice="..."` matrix in angstroms. Periodic `RG` requires the selected atoms to belong to one connected PDB bond component. See [periodic-boundary conventions](docs/periodic-boundaries.md).

## Architecture

Public headers separate trajectory/topology I/O, selections, geometry, observables, event extraction, interval algebra, AST parsing, planning, and execution. `FrameReader` is the backend-neutral streaming interface used by both the built-in XYZ reader and the optional Chemfiles XTC/TRR/DCD adapter. Observable classes are likewise independent of the parser.

The straightforward contact kernel is currently O(N×M). Its stable API permits cell lists, spatial hashing, neighbor lists, SIMD, or OpenMP underneath it. Other intended extension points are parallel frame evaluation, thread pools, memory-mapped or compressed trajectory readers, and query-plan caching. See [ROADMAP.md](ROADMAP.md).

## Current limitations

PDB supplies topology metadata and explicit `CONECT` bonds. XYZ is always supported; XTC, TRR, and DCD are available through the optional chemfiles backend. Chemfiles coordinates and cell vectors are converted from angstroms to nm, while its trajectory `time` property is already interpreted as ps. Periodic `RG` requires connected bond metadata; EnsembleQL does not yet infer standard-residue bonds, mass-weight RG, or expose a user-configurable default timestep. Event intervals use sampled timestamps and therefore do not infer behavior between frames. `AND`/`OR` combine frame predicates; temporal relations operate on extracted intervals.

## Development

The native target compiles with `-Wall -Wextra -Wpedantic`; CI additionally enables `ENSEMBLEQL_WARNINGS_AS_ERRORS`. Enable microbenchmarks with `-DENSEMBLEQL_BUILD_BENCHMARKS=ON`. The synthetic benchmarks cover contact detection at increasing atom counts, streaming event extraction, and temporal joins.

Contributions should preserve scientific definitions, add boundary-condition tests, and keep file-format backends independent from the engine.
