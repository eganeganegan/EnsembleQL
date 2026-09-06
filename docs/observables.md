# Observable definitions

## Radius of gyration

`RG(selection)` computes the unweighted root-mean-square distance from the geometric centroid:

```text
Rg² = sum_i |r_i - center|² / N
```

`RG(selection, mass_weighted=true)` instead uses [CIAAW 2024 abridged standard atomic weights](https://ciaaw.org/abridged-atomic-weights.htm) in daltons derived from PDB element symbols:

```text
center = sum_i m_i r_i / sum_i m_i
Rg² = sum_i m_i |r_i - center|² / sum_i m_i
```

The result is in nm. Every selected atom must have a known, finite, positive mass. Elements for which CIAAW does not define a standard atomic weight are rejected rather than assigned an isotope-dependent value. The default remains unweighted for backward compatibility. Both variants use the same topology-aware molecule reconstruction on periodic frames.

## Contact counts

`CONTACT_COUNT(A, B)` defaults to `mode=atom`. It counts unique unordered pairs of distinct atom indices from the two selections whose nearest-image distance is at or below the cutoff. Reverse pairs caused by overlapping selections are counted once.

`CONTACT_COUNT(A, B, mode=residue)` maps contacting atoms to residues identified by `(chain, resid)`, excludes contacts within the same residue, and counts each unique unordered residue pair once. A residue pair contributes one regardless of how many of its atom pairs are in contact.

Both modes use the 0.45 nm default cutoff and accept options in either order:

```text
CONTACT_COUNT(protein, resname ATP, mode=residue, cutoff=0.4nm)
CONTACT_COUNT(protein, resname ATP, cutoff=0.4nm, mode=atom)
```

Cutoff comparisons are inclusive and both modes honor orthorhombic and triclinic periodic cells.
