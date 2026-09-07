# Peptide–surface adsorption

This deterministic four-atom trajectory demonstrates a composable adsorption mechanism rather than presenting a physical simulation. It moves a two-site peptide through approach, contact, surface reorientation, and a stable one-picosecond adsorbed interval. The query combines PBC-aware minimum surface distance, molecular orientation, temporal succession, inclusive windows, and a duration filter.

Run it from the repository root after installing EnsembleQL:

```bash
python examples/peptide_surface_adsorption/run.py
```

The single result spans 1–4 ps. Real slab analyses should replace `resname SUR` with a topology-specific surface selection and choose cutoffs from the material, solvent, and sampling interval rather than copying the illustrative values.
