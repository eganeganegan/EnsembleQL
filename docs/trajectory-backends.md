# Trajectory backends

All trajectory decoders implement EnsembleQL's `FrameReader` interface. The query engine only receives normalized `Frame` values and does not depend on file-format libraries.

## Built-in XYZ

XYZ support has no external dependency. Coordinates and extended-XYZ lattice values are interpreted as angstroms and normalized to nm. Comment-line times are normalized to ps. Frames without a time value use the configurable fallback timestep.

## Built-in PDB

PDB support has no external dependency. A conventional PDB without `MODEL` records produces one frame. In a multi-model ensemble, each `MODEL`/`ENDMDL` block produces one frame and topology metadata is read from the first model only. Model numbers are identifiers, not physical timestamps, so frames use the configured fallback timestep.

Coordinates and `CRYST1` lengths are converted from angstroms to nm. `CRYST1` lengths and angles are retained as an orthorhombic or triclinic cell. The wwPDB unitary placeholder (`a=b=c=1 Å`, 90° angles) used for non-crystallographic structures is treated as missing cell metadata rather than a physical periodic box. Atom count, chemical identity, and ordering must remain identical across models; changes raise an error instead of silently associating coordinates with the wrong topology atoms. Atom serial numbers may continue across models and are not used as identity. Every `MODEL` must have a matching `ENDMDL`.

## Optional chemfiles

[Chemfiles 0.10](https://chemfiles.org/chemfiles/0.10.4/formats.html) or newer adds streaming readers for these coordinate-bearing formats:

- Simulation formats: Amber NetCDF (`.nc`), Amber Restart (`.ncrst`), DCD (`.dcd`), GRO (`.gro`), LAMMPS trajectories (`.lammpstrj`), Tinker ARC (`.arc`), TNG (`.tng`), TPR (`.tpr`), TRJ (`.trj`), TRR (`.trr`), and XTC (`.xtc`).
- Structure and ensemble formats: CIF (`.cif`), CML (`.cml`), CSSR (`.cssr`), mmCIF (`.mmcif`), MMTF (`.mmtf`), MOL2 (`.mol2`), Molden (`.molden`), and SDF (`.sdf`).

Uncompressed XYZ and PDB continue to use EnsembleQL's dependency-free readers. With chemfiles enabled, `.xyz` and `.pdb` files compressed using gzip, bzip2, or xz are also accepted, as are compressed ARC, CIF, CML, CSSR, GRO, LAMMPS trajectory, mmCIF, MOL2, and SDF files.

Chemfiles coordinates and cell lengths are supplied in angstroms and converted to nm. A numeric `time` frame property is interpreted as ps. If a frame has no numeric time property, EnsembleQL uses the configured fallback step, which defaults to 1 ps.

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

An ordinary build retains XYZ and PDB/ENT support if chemfiles is absent. Attempting to open a recognized optional format then raises an error describing how to enable the backend. Call `ensembleql.chemfiles_backend_available()` to inspect the installed Python extension's backend capability, or `ensembleql.supported_trajectory_extensions()` for the base-extension allowlist exposed by the build.
