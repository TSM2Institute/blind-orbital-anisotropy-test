# BOAT v2.0 — Project Initialiser

**Paste this into your new Replit agent chat as the first message.**

---

## YOUR ROLE

You are the builder for BOAT v2.0 (Blind Orbital Anisotropy Test, version 2.0) — a pre-registered blind statistical test of galaxy redshift anisotropy.

You will receive build instructions one module at a time from Graham, who orchestrates between you, Claude Opus (chief architect), and Grok (independent reviewer). You build. They design and verify. Graham conducts.

**Rules:**
- Build exactly what is specified. Do not add features, optimisations, or "improvements" unless asked.
- If something in the instructions seems wrong or ambiguous, stop and ask before building.
- All imports go at module top level, never inside functions.
- Write clean, documented Python. No clever tricks. Readability over brevity.
- Every module must be independently testable.

---

## PROJECT CONTEXT

This Replit is a clone of the BOAT v1.1 project. v1.1 has been completed (result: overall FAIL, 1 of 5 surveys passed). We are now building v2.0 with improved methodology.

**You do not need to understand the underlying cosmological theory.** Your job is to implement the blind test protocol exactly as specified. The test is honest — we publish results regardless of outcome.

---

## STEP 1: CLEAN THE CLONE

This clone contains v1.1 code and results that must be removed or replaced. Execute the following cleanup:

### DELETE these files (v1.1 artifacts, not used in v2):
- `boat/runner.py` (will be rebuilt)
- `boat/analysis.py` (will be rebuilt)
- `boat/manifest.py` (will be rebuilt)
- `boat/datasets.py` (will be rebuilt)
- `boat/results/` — delete the entire directory and its contents
- `BOAT_MANIFEST_v1.1.md` (or any v1.x manifest file)
- Any `boat_results_*.json` files in the project root
- `README.md` (will be rebuilt)
- `main.py` (will be rebuilt)

### KEEP these files (unchanged for v2):
- `boat/__init__.py`
- `boat/geometry.py` — coordinate transforms and SR Doppler (no changes needed)
- `boat/hashing.py` — SHA-256 integrity functions (no changes needed)
- `attached_assets/` — entire directory (contains survey data files)

### VERIFY the following data files exist in `attached_assets/`:
- `6dFGSzDR3.txt.gz` — 6dFGS DR3 survey data
- `MyTable_1.csv` — SDSS DR18 CasJobs export (~88 MB)

The following files may also be present but are NOT used in v2. Leave them in place, do not delete:
- `best.observations.idz.gz` (2dFGRS)
- `SpecObj.fits` (GAMA DR4)
- `asu_1771586798714.tsv` (2MRS)

### CREATE the new directory structure:
```
boat-blind-orbital-anisotropy-test-v2/
├── boat/
│   ├── __init__.py          # KEEP (existing)
│   ├── geometry.py          # KEEP (existing, unchanged)
│   ├── hashing.py           # KEEP (existing, unchanged)
│   ├── manifest.py          # TO BE BUILT (next prompt)
│   ├── montecarlo.py        # TO BE BUILT
│   ├── datasets.py          # TO BE BUILT
│   ├── analysis.py          # TO BE BUILT
│   ├── runner.py            # TO BE BUILT
│   └── results/             # CREATE empty directory
│       └── README.md        # CREATE: "Results directory. Populated by runner.py."
├── attached_assets/         # KEEP (existing, survey data)
├── main.py                  # TO BE BUILT
└── README.md                # TO BE BUILT
```

### CREATE a placeholder `README.md` in the project root:
```markdown
# BOAT v2.0 — Blind Orbital Anisotropy Test

**Status:** Under construction

**Institution:** TSM2 Institute for Cosmology Ltd

Modules are being built sequentially. Do not run until all modules are complete and the manifest is sealed.
```

### CREATE a placeholder `main.py`:
```python
"""
BOAT v2.0 — Blind Orbital Anisotropy Test
Main entry point. Do not run until all modules are built and manifest is sealed.
"""

def main():
    print("BOAT v2.0 — modules under construction")
    print("Do not execute until manifest is sealed.")

if __name__ == "__main__":
    main()
```

---

## STEP 2: VERIFY KEPT MODULES

After cleanup, run these quick checks:

### Check geometry.py imports cleanly:
```python
from boat.geometry import radec_to_unit_vector
v = radec_to_unit_vector(167.7866, -7.1454)
print(f"CMB dipole unit vector: {v}")
# Should be approximately [-0.9698, +0.2099, -0.1244]
```

### Check hashing.py imports cleanly:
```python
from boat.hashing import hash_file, hash_dataframe
import pandas as pd
df = pd.DataFrame({'a': [1, 2, 3]})
h = hash_dataframe(df)
print(f"DataFrame hash: {h[:16]}...")
# Should produce a valid hex string (first 16 chars displayed)
```

If either check fails, report the error. Do not attempt to fix geometry.py or hashing.py without instruction.

---

## STEP 3: REPORT BACK

Once cleanup and verification are complete, report:

1. Files deleted (list them)
2. Files kept (list them)
3. Data files confirmed present in attached_assets/
4. geometry.py import check: PASS/FAIL
5. hashing.py import check: PASS/FAIL
6. New directory structure confirmed

**Do not proceed to building any new modules.** Wait for the next build prompt from Graham.

---

## CRITICAL ERRORS TO AVOID (from v1.1 lessons)

These were caught during v1 development. Keep them in mind for all future modules:

1. **Residual formula:** Δz = z_obs − mean(z_obs). NEVER subtract z_model. z_model is for comparison only.
2. **Decoy seed:** Same seed for all datasets. Never SEED + i.
3. **Duplicate removal:** Use astropy `match_coordinates_sky`, not RA-sort.
4. **No threshold tuning from data.** Ratio threshold comes from Monte Carlo only.
5. **R_weight = 5° is fixed.** No calibration on real data. Geoffrey's direct instruction.
6. **z_model never enters the statistical pipeline.** It appears in output reports only.
7. **All imports at module top level.** Not inside functions.

---

*Prompt prepared by Claude Opus (Anthropic) for BOAT v2.0 build sequence.*
*21 February 2026*
