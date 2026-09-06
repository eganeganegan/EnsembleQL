# Periodic-boundary conventions

EnsembleQL currently supports orthorhombic periodic boxes for pairwise distance observables. A frame stores three positive box lengths `(Lx, Ly, Lz)` in nm. For each Cartesian displacement component, the minimum-image displacement is:

```text
delta_i = (a_i - b_i) - L_i * round((a_i - b_i) / L_i)
```

The reported distance is the Euclidean norm of that displacement. `DISTANCE`, `CONTACT`, and `CONTACT_COUNT` apply this convention automatically whenever a frame declares a box. Coordinates do not need to be wrapped into the primary cell.

## XYZ metadata

Two encodings are accepted:

```text
time=10ps box=20,20,30A
```

The compact form requires one explicit unit (`A`, `angstrom`, or `nm`) shared by all three lengths.

Extended XYZ diagonal lattice matrices are interpreted in angstroms:

```text
Lattice="20 0 0 0 20 0 0 0 30"
```

Nonzero off-diagonal components identify a triclinic cell and are rejected. Missing, malformed, non-finite, zero, or negative box lengths are also rejected.

## Whole-molecule quantities

Radius of gyration cannot generally be recovered by applying pairwise minimum images around a naive coordinate centroid. It requires bonded topology or an explicit molecule-unwrapping policy. Therefore, `RG` on a periodic frame raises an error instead of returning a potentially incorrect value.

Triclinic minimum images and topology-aware molecule reconstruction remain roadmap items. Any future implementation must retain explicit box metadata and document the chosen wrapping and unwrapping conventions.
