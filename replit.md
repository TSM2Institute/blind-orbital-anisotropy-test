# BOAT v2.0 — Blind Orbital Anisotropy Test

## Overview

BOAT v2.0 (Blind Orbital Anisotropy Test, version 2.0) is a pre-registered blind statistical test of galaxy redshift anisotropy, developed by the TSM2 Institute for Cosmology Ltd. The project tests whether observed galaxy redshifts show a dipole anisotropy aligned with a predicted axis (the CMB dipole direction), compared against randomly generated decoy axes.

**Current state:** The repository is a clone of the completed BOAT v1.1 project (which yielded an overall FAIL result — 1 of 5 surveys passed). It is being rebuilt as v2.0 with improved methodology. The project is under construction — modules are being built sequentially according to instructions from an external architect (Claude Opus) and reviewer (Grok). The user (Graham) orchestrates between the agents.

**Key principle:** This is an honest, blind scientific test. Results are published regardless of outcome. The test protocol is sealed before execution, and all parameters are locked in a manifest file.

## User Preferences

Preferred communication style: Simple, everyday language.

Additional working rules from the project owner:
- Build exactly what is specified. Do not add features, optimizations, or "improvements" unless asked.
- If something in the instructions seems wrong or ambiguous, stop and ask before building.
- All imports go at module top level, never inside functions.
- Write clean, documented Python. No clever tricks. Readability over brevity.
- Every module must be independently testable.
- Do not run the test until all modules are complete and the manifest is sealed.

## System Architecture

### Project Structure

The project is a Python package (`boat/`) with a modular design. Each module has a single responsibility:

```
boat/
  __init__.py        — Package init
  manifest.py        — Locked constants (single source of truth, no magic numbers elsewhere)
  geometry.py        — Coordinate transforms and SR Doppler predictions (no observed data)
  datasets.py        — Survey data loading, filtering, quality cuts
  analysis.py        — Residual computation, dipole regression, decoy axis comparison
  hashing.py         — SHA-256 integrity verification for datasets and outputs
  runner.py          — Main test orchestrator, produces final PASS/FAIL verdict
  results/           — Output directory for test results (JSON artifacts)
main.py              — Entry point
attached_assets/     — Survey data files (large, not in git) and build prompts
```

### Module Dependency Chain

`manifest.py` → `geometry.py` → `datasets.py` → `analysis.py` → `runner.py`

- **manifest.py** is the foundation. Every physical constant, threshold, seed, and test parameter lives here. Other modules import from it. No values may be modified after the seal date.
- **geometry.py** converts RA/Dec sky positions to unit vectors, computes line-of-sight velocities, and generates model redshift predictions (`z_model`) using special relativity Doppler formula. It never touches observed redshift data.
- **datasets.py** loads galaxy survey catalogues from various file formats (gzipped ASCII, FITS, CSV, TSV), applies quality cuts (redshift range, galactic plane masking, minimum sample size, duplicate removal), and provides synthetic data generators for validation.
- **analysis.py** computes redshift residuals (observed minus mean, isolating directional signal), runs dipole/quadrupole regression against the true predicted axis and N decoy axes, ranks axes by correlation strength, and returns a binary PASS/FAIL verdict.
- **hashing.py** computes SHA-256 hashes for files, DataFrames, and dictionaries to ensure immutability and reproducibility of all test artifacts.
- **runner.py** orchestrates the full test across all selected datasets, collects per-dataset results, applies the overall pass criterion (K of N datasets must pass), and saves the final results JSON.

### Test Logic

1. For each galaxy survey dataset:
   - Load raw catalogue and apply quality cuts
   - Compute residuals: `delta_z = z_obs - mean(z_obs)` (removes isotropic Hubble flow)
   - Compute `cos(θ)` between each galaxy direction and the true test axis (CMB dipole direction)
   - Fit dipole correlation: Pearson r between `cos(θ)` and `delta_z`
   - Generate N decoy axes (deterministic from seed) and fit dipole for each
   - PASS if: true axis has `|r| > threshold`, `p < p_threshold`, AND true axis ranks #1 among all axes

2. Overall verdict: PASS if K or more of N datasets pass (defined in manifest)

### Key Design Decisions

- **Blind test protocol:** All parameters are locked before seeing real data results. The manifest is cryptographically hashed (SHA-256) for verification.
- **Decoy axes for calibration:** Random comparison axes are generated from a seeded RNG to ensure reproducibility. The true axis must rank #1 to pass — this prevents false positives from survey geometry artifacts.
- **z_model is NOT subtracted from data:** The model prediction is computed for reference only. Residuals are simply `z_obs - mean(z_obs)`. This is a deliberate methodological choice — subtracting z_model would remove the very signal being tested.
- **Deterministic reproducibility:** All random processes use fixed seeds from the manifest (dataset selection seed, decoy axes seed).

### v2.0 Build Status (22 Feb 2026)

| Module | Status |
|--------|--------|
| `boat/manifest.py` | BUILT — 8/8 verification checks passing |
| `boat/geometry.py` | UNCHANGED from v1.1 — working |
| `boat/hashing.py` | UNCHANGED from v1.1 — working |
| `boat/montecarlo.py` | BUILT — 6dFGS + SDSS calibration complete (10,000 realisations each) |
| `boat/analysis.py` | BUILT — verified on 6dFGS (1.9s runtime) |
| `boat/datasets.py` | TO BE BUILT |
| `boat/runner.py` | TO BE BUILT |
| `main.py` | placeholder, TO BE REBUILT |

Monte Carlo calibration results (provisional, 2 of 4 surveys):
- 6dFGS: p99 = 1.5776, mean ratio = 0.8630 (119,906 galaxies)
- SDSS: p99 = 1.4984, mean ratio = 0.8000 (1,575,204 galaxies)
- Combined provisional R_NULL_99 = 1.5627
- Saved to: `boat/results/montecarlo_null_calibration_partial.json`

6dFGS analysis.py verification (R_NULL_99=None, so verdict=FAIL as expected):
- True axis |r| = 0.061962, p = 3.31e-78, ratio = 1.5489, rank = #6/51
- Orbital discriminator: DETECTED on k_hat

### Critical Errors to Avoid (from v1.1 lessons)

1. Residual formula: Δz = z_obs − mean(z_obs). NEVER subtract z_model.
2. Decoy seed: Same seed for all datasets. Never SEED + i.
3. Duplicate removal: Use astropy `match_coordinates_sky`, not RA-sort.
4. No threshold tuning from data. Ratio threshold comes from Monte Carlo only.
5. R_weight = 5° is fixed. No calibration on real data.
6. z_model never enters the statistical pipeline. It appears in output reports only.
7. All imports at module top level. Not inside functions.

## External Dependencies

### Python Packages (requirements.txt)

| Package | Purpose |
|---------|---------|
| `numpy` | Array math, vector operations, coordinate transforms |
| `pandas` | DataFrame handling for survey catalogues |
| `scipy` | Statistical tests (Pearson correlation via `scipy.stats`) |
| `astropy` | FITS file reading, sky coordinate handling (`SkyCoord`), unit conversions |
| `matplotlib` | Plotting (diagnostic visualizations) |

### Survey Data Files

Large data files (~241 MB total) are **not stored in git**. They must be downloaded separately and placed in `attached_assets/`. The datasets module expects specific filenames:

| Survey | Filename | Format | Size |
|--------|----------|--------|------|
| 6dF Galaxy Survey DR3 | `6dFGSzDR3.txt_*.gz` | Gzipped ASCII | 5.3 MB |
| 2dF Galaxy Redshift Survey | `best.observations.idz_*.gz` | Gzipped ASCII | 14 MB |
| GAMA DR4 SpecObj | `SpecObjv27_*.fits` | FITS binary table | 120 MB |
| 2MASS Redshift Survey | `asu_*.tsv` | TSV (VizieR export) | 13 MB |
| SDSS DR18 | `MyTable_TSM2_0_*.csv` | CSV (CasJobs export) | 88 MB |

Download sources are documented in `attached_assets/README.md`. No external APIs are called at runtime — all data is loaded from local files.

### No Database

This project does not use a database. All data is loaded from flat files and processed in-memory using pandas DataFrames. Results are saved as JSON files.

### No Web Server / API

This is a command-line scientific computation tool. There is no web interface, no API endpoints, and no server component.