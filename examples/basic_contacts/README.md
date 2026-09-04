# Basic contact query

Use any matching PDB/XYZ pair with:

```text
FIND CONTACT(resid 17, resid 42) FOR >= 2ns;
```

XYZ coordinates are interpreted in angstroms and comment lines may contain `time=5ns` (or `fs`, `ps`, `us`).
