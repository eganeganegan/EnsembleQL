# EnsembleQL research roadmap

## Near term

- Add PDB-coordinate frames and backend adapters for XTC, TRR, and DCD, evaluating chemfiles without coupling public reader interfaces to it.
- Add periodic boxes and explicit minimum-image conventions.
- Add mass-weighted radius of gyration, residue-residue contact modes, plan explain output, and query-plan caching.
- Optimize contacts behind the existing API with neighbor lists, cell lists/spatial hashing, SIMD, and optional OpenMP/thread-pool execution.

Completed foundation work includes precedence-aware boolean selections, configurable dimensional contact cutoffs, strict sampled-time semantics, per-frame observable sharing, and Linux/macOS CI.

## Molecular observables

Hydrogen bonds, backbone dihedrals, DSSP-like secondary structure, SASA, RMSD, coordination number, salt bridges, aromatic stacking, peptide-surface distance, and molecular orientation.

Every observable must document its mathematical definition, units, selection behavior, inequality/cutoff convention, and treatment of periodic boundaries before release.

## Temporal analysis

Add `DURING`, `UNTIL`, `REPEATS`, `PRECEDES`, and `IMMEDIATELY_FOLLOWED_BY`, followed by recurrence statistics, conditional probabilities, transition matrices, motifs, recurring subsequences, temporal clustering, event graphs, and state-transition networks.

## Research applications

For IDRs: transient salt bridges, interaction switching, long-range contact formation, collapse/expansion, secondary-structure nucleation, recurrent motifs, and residue interaction lifetimes.

For peptide/material systems: represent adsorption as composable events such as:

```text
APPROACH -> CONTACT -> SURFACE_REORIENTATION -> STABLE_ADSORPTION
```

This requires surface selections, distance/orientation observables, PBC-aware slabs, and definitions that distinguish transient contact from stable adsorption.
