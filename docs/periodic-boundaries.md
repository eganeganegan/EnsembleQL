# Periodic-boundary conventions

EnsembleQL stores a periodic cell as three Cartesian lattice vectors `a`, `b`, and `c` in nm. Orthorhombic boxes are represented by diagonal vectors; triclinic cells retain their off-diagonal components.

For a Cartesian displacement `d`, EnsembleQL converts `d` to fractional coordinates using the inverse cell matrix. It rounds to a nearby lattice translation, evaluates that translation and its 26 immediate neighbors, and selects the Cartesian image with the smallest Euclidean norm. `DISTANCE`, `CONTACT`, and `CONTACT_COUNT` apply this nearest-image search automatically whenever a frame declares a cell. Coordinates do not need to be wrapped into the primary cell.

## XYZ metadata

Two encodings are accepted:

```text
time=10ps box=20,20,30A
```

The compact form requires one explicit unit (`A`, `angstrom`, or `nm`) shared by all three lengths.

Full extended XYZ lattice matrices are interpreted as the three cell vectors in angstroms:

```text
Lattice="20 0 0 0 20 0 0 0 30"
```

```text
Lattice="20 1 0 0 20 0 0 0 30"
```

Off-diagonal components identify a triclinic cell. Missing, malformed, non-finite, or singular cell matrices are rejected. Compact `box=` metadata must contain three finite positive orthorhombic lengths.

## Whole-molecule quantities

Radius of gyration cannot generally be recovered by applying pairwise minimum images around a naive coordinate centroid. On a periodic frame, EnsembleQL starts at the first selected atom and traverses the full topology bond graph. Each bonded neighbor is placed using its nearest-image displacement from the already reconstructed atom. Paths may pass through atoms outside the final selection, which permits selections such as backbone atoms from a connected molecule.

PDB topology bonds come from `CONECT` records. EnsembleQL requires all atoms included in one `RG` selection to be reachable from one bond component; otherwise it raises an error rather than choosing an arbitrary relative image for disconnected molecules. Standard-residue bond inference is not yet implemented.
