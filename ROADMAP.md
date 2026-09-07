# EnsembleQL research roadmap

## Completed performance foundation

The contact kernel provides spatial hashing for Cartesian and fractional triclinic cells, plus reusable Verlet-style neighbor lists for multi-frame observables. OpenMP builds use explicit SIMD filtering for Cartesian candidate distances, boolean searches retain scalar early exit, and SASA uses optional threaded execution over selected atoms. Ordinary builds retain a portable serial fallback.

Completed foundation work includes precedence-aware boolean selections, configurable dimensional contact cutoffs, atom- and residue-level contact counting, spatial-hash acceleration with boolean early exit for non-periodic, orthorhombic, and triclinic contacts, mass-weighted radius of gyration, strict sampled-time semantics with configurable fallbacks, topology-aware periodic molecule unwrapping, dependency-free XYZ and PDB trajectory readers, an optional backend-neutral chemfiles adapter spanning its coordinate-bearing format set, per-frame observable sharing, topology-scoped query-plan caching, plan explanation output, and Linux/macOS CI.

## Molecular observables

Implemented observables include hydrogen bonds, backbone dihedrals, backbone-geometry helix classification, Shrake-Rupley SASA, aligned or raw first-frame RMSD, coordination number, salt bridges, parallel aromatic stacking, peptide-surface distance, and molecular orientation.

Mathematical definitions, units, selection behavior, inequality/cutoff conventions, and periodic-boundary treatment are recorded in `docs/observables.md`.

## Temporal analysis

Implemented temporal analysis includes `DURING`, `UNTIL`, `REPEATS`, `PRECEDES`, and `IMMEDIATELY_FOLLOWED_BY`, plus recurrence statistics, observation-time-normalized event frequencies, conditional probabilities, transition matrices, motifs, recurring subsequences, temporal clustering, event graphs, and state-transition networks.

## Research applications

For IDRs: transient salt bridges, interaction switching, long-range contact formation, collapse/expansion, secondary-structure nucleation, recurrent motifs, and residue interaction lifetimes.

For peptide/material systems: represent adsorption as composable events such as:

```text
APPROACH -> CONTACT -> SURFACE_REORIENTATION -> STABLE_ADSORPTION
```

The implemented surface selections, distance/orientation observables, PBC-aware slab handling, and explicit duration predicates distinguish transient contact from stable adsorption.

Auditable executable fixtures now cover both application families: `examples/idr_contact_switching` exercises switching, recurrence, and event-graph analysis, while `examples/peptide_surface_adsorption` composes approach, contact, reorientation, and duration-qualified adsorption under PBC. Their synthetic cutoffs are semantic examples, not scientific defaults.

## Validation and release milestones

1. **Completed:** validate every observable family against independent reference implementations on published, versioned trajectories and record thresholds, tolerances, checksums, and provenance in machine-readable manifests.
2. Ship reproducible Linux and macOS wheels and release-level compatibility tests across supported trajectory formats.
3. Profile production IDR and slab workloads, then add bounded-memory parallel frame scheduling or accelerator kernels only where measurements show a durable gain.
