# Scientific validation

This directory contains reproducible comparisons against MDAnalysis, MDTraj, and NumPy. Generated trajectory data live under `validation/data/` and are intentionally excluded from Git. The committed JSON manifests record dataset provenance, source and prepared-file SHA-256 digests, dependency versions, exact queries and thresholds, event intervals, and time tolerances.

## Run the complete suite

Install the pinned validation dependencies and prepare both fixtures:

```bash
python -m pip install -e ".[validation,test]"
python validation/prepare_adk.py
python validation/prepare_membrane_peptide.py --download
```

Regenerate the independent reference manifests after an intentional scientific-method change, or use `--check` to prove they are current:

```bash
python validation/generate_adk_reference.py
python validation/generate_membrane_peptide_reference.py
python validation/generate_adk_reference.py --check
python validation/generate_membrane_peptide_reference.py --check
```

Run the validation tests:

```bash
pytest tests/python/test_scientific_validation.py
```

The test skips when external files or the Chemfiles backend are unavailable. It fails, rather than silently using different data, when a prepared file does not match its recorded digest. A dedicated CI job performs preparation, stale-manifest checks, and all comparisons on every push and pull request.

## Adenylate kinase

The 98-frame, 3,341-atom adenylate kinase DIMS trajectory is distributed with MDAnalysisTests 2.10.0. It is associated with the closed-to-open transition reported in [Beckstein et al., *Journal of Molecular Biology* 394 (2009), 160–176](https://doi.org/10.1016/j.jmb.2009.09.009).

Preparation writes a PDB topology with explicit `CONECT` bonds and converts the original DCD to XTC. The XTC is used for scientific comparisons because its timestamps agree with MDAnalysis; the original DCD is retained for format-compatibility testing.

The AdK manifest validates `DISTANCE`, `CONTACT`, atom- and residue-mode `CONTACT_COUNT`, unweighted and mass-weighted `RG`, `HBOND`, `DIHEDRAL`, `HELIX`, `SASA`, aligned and raw `RMSD`, `COORDINATION_NUMBER`, `SALT_BRIDGE`, and `AROMATIC_STACKING`.

## Membrane peptide

The periodic fixture is derived from the published [MDAnalysis membrane/peptide tutorial dataset](https://doi.org/10.6084/m9.figshare.8046437), licensed CC BY 4.0. It retains the 197-atom peptide, the 128 DMPC phosphorus atoms, every tenth frame, and each frame's changing orthorhombic unit cell. The resulting fixture has 325 atoms and 101 frames at 100 ps spacing.

The membrane manifest validates periodic `SURFACE_DISTANCE` against MDAnalysis pair distances and periodic `ORIENTATION` against MDAnalysis minimum-image vectors. Together, the two manifests cover every currently implemented observable family.

## Interpretation

These are event-level conformance checks: independent software computes a Boolean state at every sampled frame using the recorded threshold, then EnsembleQL must return the identical maximal time intervals and matching-frame count. The time tolerance applies only to serialized trajectory timestamps. The membrane manifest also records the smallest numerical distance from any reference value to its decision threshold so boundary sensitivity is auditable.
