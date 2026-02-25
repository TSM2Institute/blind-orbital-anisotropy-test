# BOAT v2.0 — Blind Orbital Anisotropy Test

**Status:** Production run complete — overall verdict: **FAIL**

**Institution:** TSM2 Institute for Cosmology Ltd  
**Repository:** [github.com/TSM2Institute/blind-orbital-anisotropy-test](https://github.com/TSM2Institute/blind-orbital-anisotropy-test)  
**Predecessor:** BOAT v1.1 (sealed 19 February 2026, overall verdict FAIL)

---

## What Is BOAT?

BOAT (Blind Orbital Anisotropy Test) is a pre-registered blind statistical test that checks whether galaxy redshift residuals show a dipole anisotropy aligned with the CMB dipole velocity direction, as predicted by the TSM2.1 orbital motion model.

All parameters, thresholds, and datasets are locked before any real data touches the analysis pipeline. Results are published regardless of outcome.

## v2.0 Results (25 February 2026)

**Manifest sealed:** 23 February 2026  
**R_NULL_99 = 1.541144** (Monte Carlo calibrated from 40,000 null realisations across 4 surveys)

| Survey | N galaxies | \|r\| | Ratio | Rank | p-value | Ratio Pass | Rank Pass | z-bin Pass | Verdict |
|--------|-----------|-------|-------|------|---------|------------|-----------|------------|---------|
| 6dFGS DR3 | 119,906 | 0.0620 | 1.5489 | #6/51 | 3.31e-78 | PASS | FAIL | FAIL | **FAIL** |
| SDSS DR18 | 1,575,204 | 0.0827 | 1.4655 | #5/51 | 2.03e-05 | FAIL | FAIL | FAIL | **FAIL** |
| DESI DR1 BGS | 5,841,593 | 0.0020 | 0.5313 | #36/51 | 2.10e-04 | FAIL | FAIL | FAIL | **FAIL** |
| DESI DR1 LRG | 3,772,804 | 0.0223 | 1.0616 | #27/51 | 1.12e-205 | FAIL | FAIL | PASS | **FAIL** |

**Overall: FAIL** (0 of 4 passed, 3 required)

### Orbital Discriminator (Secondary Test)

After subtracting the primary CMB dipole, residuals were tested against the CoG X direction (r̂) and the orbital normal (k̂):

| Survey | Axis | Post-subtraction p | Ratio | Detection |
|--------|------|-------------------|-------|-----------|
| 6dFGS | k̂ | 4.02e-07 | 1.853 | **DETECTED** |
| SDSS | — | — | — | NOT DETECTED |
| DESI BGS | r̂ | 2.42e-31 | 1.885 | **DETECTED** |
| DESI LRG | r̂ | 0.00e+00 | 2.334 | **DETECTED** |

**Orbital signal detected in 3 of 4 surveys.**

## v1.1 Results (21 February 2026)

Overall verdict: **FAIL** (1 of 5 passed, 3 required). Only 6dFGS passed. See v1.1 tag for details.

## Methodology

1. Compute residuals: Δz = z_obs − mean(z_obs)
2. Fit dipole correlation between Δz and cos(θ) to the CMB dipole axis
3. Compare against 50 seeded random decoy axes
4. PASS requires: ratio ≥ R_NULL_99, rank in top 3, p < 0.01, z-bin invariance
5. Overall PASS: ⌈N/2⌉ + 1 of N surveys must pass

See `BOAT_v2_MANIFEST.md` for the complete pre-registered protocol.

## Artifacts

| File | Description |
|------|-------------|
| `BOAT_v2_MANIFEST.md` | Sealed manifest (SHA-256: e99156a9...) |
| `boat/results/boat_results_production_v2.0.json` | Production results (SHA-256: c841dc85...) |
| `boat/results/montecarlo_null_calibration_final.json` | MC calibration (SHA-256: d661f933...) |

## Project Structure

```
boat/
  manifest.py      — Locked constants (R_NULL_99 = 1.541144)
  geometry.py      — Coordinate transforms, SR Doppler predictions
  datasets.py      — Survey loading, quality cuts, pre-qualification
  analysis.py      — Dipole regression, decoy comparison, z-bin test
  montecarlo.py    — Null calibration (40,000 realisations)
  runner.py        — Production pipeline orchestrator
  hashing.py       — SHA-256 integrity verification
  results/         — Output JSON artifacts
```

## Authors

**Theory:** Geoffrey E. Thwaites — TSM2.1 Author  
**Implementation & Coordination:** Graham Hill — TSM2 Institute  
**Verification:** Claude (Anthropic), Grok (xAI)

---

*"The simplest test is the most decisive."*

*Pipeline is honest. We record what we see.*
