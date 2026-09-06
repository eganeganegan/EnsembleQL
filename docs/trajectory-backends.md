# Trajectory backends

All trajectory decoders implement EnsembleQL's `FrameReader` interface. The query engine only receives normalized `Frame` values and does not depend on file-format libraries.

## Built-in XYZ

XYZ support has no external dependency. Coordinates and extended-XYZ lattice values are interpreted as angstroms and normalized to nm. Comment-line times are normalized to ps.

## Optional chemfiles

Chemfiles 0.10 or newer adds streaming readers for XTC, TRR, and DCD. Its coordinates and cell lengths are supplied in angstroms and converted to nm. The `time` frame property for these formats is supplied in ps. If a file has no time property, EnsembleQL uses the reader's default 1 ps step.

The adapter accepts infinite, orthorhombic, and triclinic cells. Chemfiles cell matrices are converted into EnsembleQL's three-vector cell representation without discarding cell angles.

Configure against an installed package:

```bash
cmake -S . -B build -DENSEMBLEQL_REQUIRE_CHEMFILES=ON
```

Or fetch the pinned stable 0.10.4 source:

```bash
cmake -S . -B build \
  -DENSEMBLEQL_FETCH_CHEMFILES=ON \
  -DENSEMBLEQL_REQUIRE_CHEMFILES=ON
```

For an editable Python installation, pass the same options through scikit-build:

```bash
CMAKE_ARGS="-DENSEMBLEQL_FETCH_CHEMFILES=ON -DENSEMBLEQL_REQUIRE_CHEMFILES=ON" \
  python -m pip install -e .
```

An ordinary build remains dependency-free if chemfiles is absent. Attempting to open XTC, TRR, or DCD then raises an error describing how to enable the backend. Call `ensembleql.chemfiles_backend_available()` to inspect the installed Python extension's capability.
