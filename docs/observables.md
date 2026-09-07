# Observable definitions

## Distances and contacts

`DISTANCE(A, B)` returns the minimum distance in nm over all distinct atom pairs drawn from the two selections. `CONTACT(A, B)` is true when that value is at or below its inclusive cutoff, which defaults to 0.45 nm. Both reject selections that contain no distinct pair. Cartesian distance is used without a cell; declared orthorhombic and triclinic cells use nearest-image distance.

`CONTACT` and the contact counts below use Verlet-style reuse during trajectory queries. This changes candidate generation only, never the cutoff or reported distance.

## Radius of gyration

`RG(selection)` computes the unweighted root-mean-square distance from the geometric centroid:

$$
R_g^2 = \frac{1}{N}\sum_{i=1}^{N}
\left\lVert \mathbf{r}_i - \mathbf{r}_{\mathrm{center}} \right\rVert^2
$$

`RG(selection, mass_weighted=true)` instead uses [CIAAW 2024 abridged standard atomic weights](https://ciaaw.org/abridged-atomic-weights.htm) in daltons derived from PDB element symbols:

$$
\begin{aligned}
\mathbf{r}_{\mathrm{center}}
  &= \frac{\sum_{i=1}^{N} m_i\mathbf{r}_i}{\sum_{i=1}^{N}m_i}, \\
R_g^2
  &= \frac{\sum_{i=1}^{N}m_i
      \left\lVert\mathbf{r}_i-\mathbf{r}_{\mathrm{center}}\right\rVert^2}
      {\sum_{i=1}^{N}m_i}.
\end{aligned}
$$

The result is in nm. Every selected atom must have a known, finite, positive mass. Elements for which CIAAW does not define a standard atomic weight are rejected rather than assigned an isotope-dependent value. The default remains unweighted for backward compatibility. Both variants use the same topology-aware molecule reconstruction on periodic frames.

## Contact counts

`CONTACT_COUNT(A, B)` defaults to `mode=atom`. It counts unique unordered pairs of distinct atom indices from the two selections whose nearest-image distance is at or below the cutoff. Reverse pairs caused by overlapping selections are counted once.

`CONTACT_COUNT(A, B, mode=residue)` maps contacting atoms to residues identified by `(chain, resid)`, excludes contacts within the same residue, and counts each unique unordered residue pair once. A residue pair contributes one regardless of how many of its atom pairs are in contact.

Both modes use the 0.45 nm default cutoff and accept options in either order:

```text
CONTACT_COUNT(protein, resname ATP, mode=residue, cutoff=0.4nm)
CONTACT_COUNT(protein, resname ATP, cutoff=0.4nm, mode=atom)
```

Cutoff comparisons are inclusive and both modes honor orthorhombic and triclinic periodic cells. Spatial hashing accelerates sufficiently large searches without changing pair order or counting semantics; triclinic cells are partitioned in wrapped fractional coordinates, while small selections use the pairwise reference path. Stateful observables reuse candidates built at cutoff plus a 0.1 nm skin until an atom exceeds half-skin displacement or the cell changes. Boolean `CONTACT` uses the same cutoff convention but stops searching after the first matching atom pair.

## Hydrogen bonds

`HBOND(donors, acceptors)` is true when at least one selected donor has an explicitly bonded hydrogen and an acceptor satisfying both criteria:

$$
d(D,A) \le 0.35\ \mathrm{nm},
\qquad
\angle D\!H\!A \ge 150^{\circ}.
$$

Override these inclusive defaults with `distance=` and `min_angle=`. Hydrogens are obtained from PDB `CONECT` bonds; EnsembleQL does not guess covalent connectivity. All displacements use nearest periodic images.

## Dihedrals and helix classification

`DIHEDRAL(A, B, C, D)` requires four selections that each resolve to exactly one atom. It returns the signed torsion angle in the range $[-180^{\circ}, 180^{\circ}]$, using nearest-image bond vectors.

`HELIX(selection)` is a documented backbone-geometry classifier, not a full DSSP implementation. For every selected residue with complete neighboring $\mathrm{C}_{i-1}$–$\mathrm{N}_i$–$\mathrm{CA}_i$–$\mathrm{C}_i$–$\mathrm{N}_{i+1}$ atoms, it classifies the residue as alpha-helical when:

$$
-100^{\circ} \le \phi_i \le -30^{\circ},
\qquad
-80^{\circ} \le \psi_i \le -5^{\circ}.
$$

The observable is true when the classified fraction is at least `minimum_fraction`, which defaults to 0.5. Residues are kept chain-local and must have consecutive residue numbers.

## Solvent-accessible surface area

`SASA(selection)` implements the [Shrake-Rupley point-sampling definition](https://pubmed.ncbi.nlm.nih.gov/4760134/). For each selected atom $i$, a sphere with expanded radius $r_i=r_i^{\mathrm{vdW}}+r_{\mathrm{probe}}$ is sampled using a deterministic Fibonacci lattice. A point is accessible when it is outside every other atom's expanded sphere. The returned area is:

$$
\mathrm{SASA}
= \sum_i 4\pi r_i^2
  \frac{N_{i,\mathrm{accessible}}}{N_{\mathrm{points}}}.
$$

The default water probe is 0.14 nm and the default resolution is 96 points per atom. Use `probe=` and `points=` to change them. H, C, N, O, F, Si, P, S, Cl, Br, and I radii are built in; unknown elements are rejected during planning. All topology atoms occlude solvent, even when only a subset contributes reported area, and occlusion uses nearest-image distance when a cell is declared. Results are in nm² and comparisons accept `nm2`, `A2`, or `angstrom2`.

## RMSD

`RMSD(selection)` uses the first trajectory frame as the reference, centers both selected coordinate sets, and applies the [Kabsch least-squares proper rotation](https://doi.org/10.1107/S0567739476001873) before calculating:

$$
\mathrm{RMSD}
= \sqrt{\frac{1}{N}\sum_{i=1}^{N}
  \left\lVert \mathbf{R}\widetilde{\mathbf{r}}_i
  - \widetilde{\mathbf{r}}_i^{\,(0)} \right\rVert^2},
$$

where the tildes denote coordinates translated to their respective geometric centers, $\mathbf{R}$ is the fitted proper rotation, and superscript $(0)$ denotes the first-frame reference.

Set `align=false` to retain rigid-body translation and rotation. The result is in nm. Periodic selections are unwrapped around their first selected atom before comparison. Every RMSD observable is evaluated on every frame, including inside boolean expressions, so short-circuiting cannot shift its reference frame.

## Coordination and salt bridges

`COORDINATION_NUMBER(centers, neighbors)` returns the mean number of unique selected neighbors per selected center at or below an inclusive cutoff, defaulting to 0.35 nm. It is dimensionless and honors periodic cells.

`SALT_BRIDGE(positive, negative)` is true when the explicitly selected positive and negative charge-center atoms are within the inclusive cutoff, defaulting to 0.4 nm. Charge assignment is deliberately left to the selections; EnsembleQL does not silently infer protonation states.

## Aromatic stacking

`AROMATIC_STACKING(ringA, ringB)` implements parallel stacking. Each ring selection must contain at least three consistently ordered, non-collinear atoms. The predicate requires centroid distance at or below 0.55 nm and the unoriented angle between polygon normals at or below 30 degrees. Override these with `distance=` and `max_angle=`. Centroids, normals, and separation honor periodic cells.

## Surface distance and orientation

`SURFACE_DISTANCE(molecule, surface)` returns the minimum nearest-image atom distance in nm. It is intended for peptide/material approach and adsorption predicates without assuming a particular surface chemistry.

`ORIENTATION(origin, target, axis=z)` computes the unoriented angle in degrees between the nearest-image vector from the first selection's centroid to the second selection's centroid and Cartesian `x`, `y`, or `z`. Its range is $[0^{\circ}, 90^{\circ}]$; reversing the molecular axis does not change the result.
