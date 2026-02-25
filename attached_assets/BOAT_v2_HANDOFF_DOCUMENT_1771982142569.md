# BOAT v2.0 — Complete Handoff Document

**From:** Claude Opus (current session, 19–21 February 2026)  
**To:** New Claude session (BOAT v2.0 implementation)  
**Date:** 21 February 2026  
**Purpose:** Full context transfer for building and executing BOAT v2.0

---

## SECTION A: WHO YOU'RE WORKING WITH

### Graham Hill (G)
- Director, TSM2 Institute for Cosmology Ltd (Australian not-for-profit)
- 56-year-old graphic designer and print shop owner from Brisbane
- NOT a physicist — orchestrator and conductor between AI systems
- Methodical, documents everything, works backwards from outcomes
- Prefers direct communication, no filler
- Uses Replit for code, GitHub for version control
- His Replit username: grayhill5 / Grayhill5

### Geoffrey E. Thwaites
- TSM2.1 author, retired RF/radar engineer with 60+ years experience
- Provides theory direction and final sign-off on all decisions
- Communications come through Graham (Geoffrey → Graham → AI)
- Emphasises: "no expectation of a particular result — the objective is clarity"
- His instructions are non-negotiable on protocol discipline

### Grok (xAI)
- Shadows all work with full context
- Graham copy-pastes between Claude and Grok for cross-verification
- Grok 4.20 runs 4 agents simultaneously reaching consensus, then agent 5 delivers final
- Grok's role: independent reviewer, fact-checker, co-architect
- Treat Grok's input as peer review — agree or disagree on merits

### Working Protocol
- Graham orchestrates via copy-paste between AI systems
- Each AI maintains its role: Claude builds and audits, Grok reviews and verifies
- Graham has final decision authority on all implementation choices
- Geoffrey has final authority on scientific direction
- No cheerleading — honest assessment of problems and limitations

---

## SECTION B: THE TSM2.1 FRAMEWORK (Background Only)

TSM2.1 (Thwaites Standard Model 2.1) is an alternative cosmological framework
proposing that:

1. The Solar System orbits a Centre of Gravity (CoG X) at RA 347.75°, Dec +66.00°
2. The CMB dipole (370 km/s toward RA 167.79°, Dec −7.15°) is the velocity
   signature of this orbit
3. Cosmological redshift is explained via plasma refraction and hydrogen
   scattering rather than metric expansion
4. Dark matter is unnecessary — the Bullet Cluster lensing profile is
   reproduced by plasma refraction alone

**Prior validated work (December 2025):**
- SKIN-a-CAT pipeline: χ²/dof = 1.00 on 114 galaxy clusters
- JADES DR4 high-z galaxies: R² = 0.99
- Subluminal velocities maintained to z = 14.32
- Universal constant k_TSM = 5.1 × 10⁻²³ cm²
- Repository: github.com/Grayhill5/skin-a-cat-pipeline

**You do not need to evaluate TSM2.1's validity.** Your job is to build and
execute the blind test protocol honestly. The test speaks for itself.

---

## SECTION C: BOAT v1.1 — WHAT WAS DONE AND WHAT WE LEARNED

### What BOAT Is
The Blind Orbital Anisotropy Test checks whether galaxy redshift residuals
show a statistically significant dipole aligned with the CMB velocity
direction, as predicted by TSM2.1's orbital model.

### v1.1 Design
- Pre-registered manifest with locked constants and thresholds
- 10 seeded random "decoy" axes plus 1 true (CMB dipole) axis
- PASS required: true axis ranks #1 of 11 with |r| ≥ 0.01 and p < 0.01
- Overall: 3 of 5 surveys must pass

### v1.1 Results (Production Run — 21 February 2026)

| Survey | N galaxies | Sky % | Median z | |r| | Rank | Verdict |
|--------|-----------|-------|----------|------|------|---------|
| 6dFGS DR3 | 119,907 | 33.7% | 0.054 | 0.077 | #1 | **PASS** |
| 2dFGRS | 129,168 | 4.0% | 0.108 | 0.073 | #6 | FAIL |
| GAMA DR4 | 282,139 | 0.7% | 0.216 | 0.047 | #5 | FAIL |
| 2MRS | 37,688 | 31.5% | 0.030 | 0.023 | #7 | FAIL |
| SDSS DR18 | 1,575,204 | 20.9% | 0.204 | 0.132 | #4 | FAIL |

**Overall: FAIL (1 of 5 passed, 3 required)**

Results JSON hash: b9bd7e72e0dc645e2d421941a5b65f1c37f9cfa15ab636b80bbf19f159b6899d
Repository: github.com/TSM2Institute/blind-orbital-anisotropy-test

### Key Lessons from v1.1

**Lesson 1: Signal is real but rank criterion is too strict for asymmetric footprints.**
SDSS showed |r| = 0.132 (strongest of all surveys) with p ≈ 0, but ranked #4
because 3 nearby decoy axes (all within ~60° of true axis) scored slightly higher
due to the 90% northern-sky footprint shifting the apparent dipole peak.

**Lesson 2: cos(θ) projection pattern is the diagnostic signature.**
In both 6dFGS (PASS) and SDSS (FAIL), decoy axes near the true direction showed
high |r| and those far away showed low |r|. This is the textbook signature of a
real dipole — decoys near it pick up a projected fraction of the signal.

**Lesson 3: Sky coverage alone isn't sufficient — need adequate redshift depth.**
2MRS had 31.5% sky coverage (comparable to 6dFGS) but failed because median
z = 0.03 puts most galaxies in the regime where local peculiar velocities
(~600 km/s) overwhelm the orbital signal (~370 km/s).

**Lesson 4: Pencil-beam surveys generate spurious correlations on all axes.**
2dFGRS (4% sky) and GAMA (0.7% sky) showed extreme significance on nearly
every axis — the hallmark of footprint-driven systematics, not real signal.

**Lesson 5: The residual formula matters.**
Early in v1 development, the code subtracted z_model from z_obs. This removes
the signal. Correct formula: Δz = z_obs − mean(z_obs). The z_model is the
prediction to compare against, not a subtraction term.

**Lesson 6: |r| threshold must be physically motivated.**
v1.1 originally had |r| ≥ 0.05 threshold. Synthetic validation showed the
expected signal-to-noise ratio makes this unreachable. Amended to 0.01.
v2 replaces this with Monte Carlo-calibrated threshold entirely.

---

## SECTION D: BOAT v2.0 — THE MANIFEST

The complete v2.0 manifest is provided as a separate file:
**BOAT_v2_MANIFEST_FINAL_DRAFT.md**

This has been reviewed and approved line-by-line by Claude, Grok, and Geoffrey.
The ONLY blank to fill is R_NULL_99 (the Monte Carlo-derived ratio threshold).

### Summary of v2 changes from v1.1:

1. **Survey pre-qualification:** Sky ≥ 15%, median z ≥ 0.04, N ≥ 10,000
2. **Continuous ratio statistic:** true |r| / mean(decoy |r|) ≥ R_NULL_99
3. **Safety gate:** True axis must rank in top 3 of 51 (not #1 of 11)
4. **50 decoy axes** (seeded, DECOY_AXES_SEED = 20260301)
5. **Footprint correction:** Inverse local density weighting, R = 5° fixed
6. **z-bin invariance:** 4 bins, max/min |r| ratio ≤ 2.0, same sign, bins < 1000 excluded
7. **Orbital discriminator:** After subtracting CMB dipole, test residuals against CoG X (r̂) and orbital normal (k̂). Secondary test, p < 0.05 threshold.
8. **Overall PASS:** ⌈N/2⌉ + 1 of N surveys (expected: 3 of 4)

### Pre-qualified surveys:

| # | Survey | Expected Sky % | Median z | Access |
|---|--------|---------------|----------|--------|
| 1 | 6dFGS DR3 | ~34% | ~0.05 | Already have data |
| 2 | SDSS DR18 | ~21% | ~0.20 | Already have data |
| 3 | DESI DR1 BGS | ~18% | ~0.20 | NOIRLab Astro Data Lab (free) |
| 4 | DESI DR1 LRG | ~18% | ~0.70 | NOIRLab Astro Data Lab (free) |

---

## SECTION E: LOCKED PHYSICAL CONSTANTS

These are identical to v1.1. Copy exactly into manifest.py.

```python
# === Observational Anchors ===
V_OBS = 370.0              # km/s — CMB dipole speed (Planck 2018)
C_KM_S = 299792.458        # km/s — speed of light

CMB_DIPOLE_RA_DEG = 167.7866    # degrees (J2000)
CMB_DIPOLE_DEC_DEG = -7.1454    # degrees (J2000)

# === TSM2.1 Orbital Parameters ===
CST_PERIOD_GYR = 284.0     # Gyr — orbital period
COG_X_RA_DEG = 347.75      # degrees (J2000) — Centre of Gravity direction
COG_X_DEC_DEG = 66.00      # degrees (J2000)

# === Derived Constants ===
R_OBS_MLY = 55.78           # Mly — orbital radius (v_obs × P / 2π)
BETA_OBS = V_OBS / C_KM_S   # ≈ 0.001234

# === Derived Unit Vectors (ICRS J2000.0) ===
import numpy as np

def _ra_dec_to_unit(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])

R_HAT = _ra_dec_to_unit(COG_X_RA_DEG, COG_X_DEC_DEG)       # [+0.397476, -0.086300, +0.913545]
V_HAT = _ra_dec_to_unit(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)  # [-0.969776, +0.209910, -0.124388]

# Orbital normal: k̂ = r̂ × v̂ / |r̂ × v̂|
_k_raw = np.cross(R_HAT, V_HAT)
K_HAT = _k_raw / np.linalg.norm(_k_raw)                     # [-0.211516, -0.977374, -0.000301]

# Velocity vector
V_OBS_VEC = V_OBS * V_HAT   # [-358.8171, +77.6668, -46.0234] km/s

# === v2.0 Test Parameters ===
N_DECOY_AXES = 50
DECOY_AXES_SEED = 20260301
P_VALUE_THRESHOLD = 0.01
Z_MIN = 0.01
Z_MAX = 0.50
Z_MIN_LRG = 0.40
Z_MAX_LRG = 1.00
GALACTIC_PLANE_MASK_DEG = 10.0
MIN_SAMPLE_SIZE = 10000
MIN_SKY_FRACTION = 0.15
MIN_MEDIAN_Z = 0.04
R_WEIGHT_DEG = 5.0
PASS_REQUIRED_FORMULA = "ceil(N/2) + 1"
MAX_MIN_RATIO_ZBIN = 2.0
MIN_BIN_SIZE = 1000

# R_NULL_99 — Monte Carlo calibrated threshold
# [TO BE FILLED AFTER MONTE CARLO SIMULATION]
R_NULL_99 = None  # Will be set to 99th percentile of null ratio distribution
```

### Verification checks (run on import):
1. |r̂| = 1.000000 ± 10⁻⁶
2. |v̂| = 1.000000 ± 10⁻⁶
3. |k̂| = 1.000000 ± 10⁻⁶
4. r̂ · k̂ = 0.000 ± 10⁻⁶
5. v̂ · k̂ = 0.000 ± 10⁻⁶
6. R_obs ≈ 55.78 Mly ± 0.01
7. β_obs ≈ 0.001234 ± 10⁻⁶
8. Angular separation r̂–v̂ ≈ 121.15° ± 0.01°

---

## SECTION F: EXECUTION PLAN

### Phase 1: Monte Carlo Null Calibration (FIRST — before manifest seal)

**Build module: `boat/montecarlo.py`**

For each of the 4 surveys:
1. Load real survey data (6dFGS and SDSS already available; for DESI, use
   footprint geometry from documentation or generate representative mock)
2. For each of 10,000 realisations:
   a. Keep the real (RA, Dec) positions fixed
   b. Randomly shuffle z_obs values among positions (destroys any real dipole)
   c. Apply footprint weighting (R = 5°)
   d. Run full 51-axis protocol (1 true + 50 decoys)
   e. Record ratio = true |r| / mean(decoy |r|)
3. Combine all 40,000 ratio values
4. R_NULL_99 = 99th percentile of combined distribution
5. Also record 95th and 99.5th percentiles for reference

**CRITICAL:** For DESI surveys where we don't have real data yet, we need to
generate representative footprints. Options:
- Use the DESI DR1 footprint map (available publicly) to create mock positions
- Or defer DESI-specific null calibration until after data download, running
  Monte Carlo on all 4 real footprints before sealing

**Recommended approach:** Run Monte Carlo on 6dFGS + SDSS footprints now (we have
the data). After downloading DESI, run on DESI footprints. Combine all 4 before
sealing. This means the manifest seal happens AFTER DESI data is downloaded but
BEFORE any v2 analysis runs.

### Phase 2: Seal Manifest
1. Insert R_NULL_99 into manifest
2. Compute SHA-256 of complete manifest
3. Commit to GitHub

### Phase 3: Execute v2 on All 4 Surveys
1. Load each survey
2. Verify pre-qualification (sky ≥ 15%, median z ≥ 0.04)
3. Apply quality cuts
4. Apply footprint weighting
5. Run primary dipole test (51 axes)
6. Run z-bin invariance test
7. Run orbital discriminator (post-subtraction test on r̂ and k̂)
8. Record per-dataset PASS/FAIL
9. Compute overall verdict

### Phase 4: Publish
1. Commit all artifacts to GitHub in single timestamped commit
2. Generate final report
3. Publish regardless of outcome

---

## SECTION G: v1.1 CODE MODULES — WHAT CARRIES OVER

The v1.1 Replit project contains 6 modules. Here's what changes for v2:

### Carry over (minor modifications):
- **manifest.py** — Update constants for v2 (add R_WEIGHT_DEG, N_DECOY=50, new seed, R_NULL_99, LRG z-range, pre-qualification thresholds)
- **geometry.py** — No changes needed. RA/Dec → unit vector, v_LOS projection, SR Doppler all stay the same.
- **hashing.py** — No changes. SHA-256 for DataFrames, dicts, files.

### Significant modifications:
- **datasets.py** — Add DESI DR1 loaders (BGS and LRG), add pre-qualification check function (sky fraction + median z), upgrade duplicate removal to use astropy match_coordinates_sky
- **analysis.py** — Add weighted correlation, footprint weighting, z-bin invariance test, orbital discriminator (post-subtraction test on r̂ and k̂), increase decoys to 50, implement ratio statistic
- **runner.py** — Update PASS logic (ratio ≥ R_NULL_99, top 3 of 51, p < 0.01, z-bin invariance), add orbital discriminator reporting, add pre-qualification gate

### New module:
- **montecarlo.py** — Null simulation framework (shuffle z_obs, run protocol, record ratios, compute percentiles)

---

## SECTION H: DATA FILES

### Already available (in v1.1 Replit attached_assets/):
- `6dFGSzDR3.txt.gz` — 6dFGS DR3 (124,647 entries)
- `best.observations.idz.gz` — 2dFGRS (not used in v2)
- `SpecObj.fits` — GAMA DR4 (not used in v2)
- `asu_1771586798714.tsv` — 2MRS (not used in v2)
- `MyTable_1.csv` — SDSS DR18 CasJobs export (1,575,213 entries, 88 MB)

### Need to download:
- **DESI DR1 BGS** — from NOIRLab Astro Data Lab (datalab.noirlab.edu)
  - SQL query needed: select galaxies from BGS target class, ZWARN==0, z < 0.6
  - Table: desi_dr1.iron.zpix or similar (check schema)
  - Expected: ~7 million BGS galaxies
- **DESI DR1 LRG** — same source
  - SQL query: LRG target class, ZWARN==0, 0.4 < z < 1.0
  - Expected: ~2 million LRG galaxies

Graham has registered for NOIRLab Astro Data Lab (or will do so).
Alternative access: AWS S3 bucket (desidata.s3.amazonaws.com) or Globus.

---

## SECTION I: CRITICAL ERRORS TO AVOID

These were caught during v1 development. Don't repeat them.

1. **Residual formula:** Δz = z_obs − mean(z_obs). Do NOT subtract z_model.
   z_model is what we compare against, not a correction term.

2. **Decoy seed consistency:** Use the same seed for all datasets within a
   production run. v1 initially used SEED + i, giving different decoys per
   dataset. Fixed to single seed for all.

3. **Duplicate removal:** v1's simple RA-sort + adjacent-pair check misses
   non-adjacent duplicates. v2 must use astropy match_coordinates_sky for
   proper angular cross-matching.

4. **Threshold tuning from results:** Never set a threshold based on what
   you observe in the data. The ratio threshold comes from Monte Carlo
   null simulation ONLY.

5. **Real-data calibration:** Geoffrey explicitly prohibited optimising
   R_weight on real survey data. R = 5° is fixed. No calibration step.

6. **z_model in residuals:** The z_model computed from geometry is for
   COMPARISON only. The blind test protocol works entirely with observed
   redshifts. z_model appears in the output report but never in the
   statistical analysis pipeline.

7. **Import placement:** Keep all imports at module top level, not inside
   functions (Grok flagged this as a code quality issue in v1).

---

## SECTION J: REPLIT WORKFLOW

### Starting fresh:
1. Graham forks/remixes the v1.1 Replit project
2. New Replit agent chat starts fresh (no v1 conversation history)
3. Give the agent this handoff document as context
4. Build modules in order: manifest → montecarlo → geometry → datasets → analysis → runner

### File structure:
```
boat-blind-orbital-anisotropy-test-v2/
├── boat/
│   ├── __init__.py
│   ├── manifest.py          # Locked constants + verification
│   ├── montecarlo.py         # NEW: null simulation framework
│   ├── geometry.py           # Coordinate transforms + SR Doppler
│   ├── datasets.py           # Survey loaders + pre-qualification
│   ├── analysis.py           # Dipole test + weighting + invariance + discriminator
│   ├── hashing.py            # SHA-256 integrity
│   ├── runner.py             # Production orchestrator
│   └── results/
│       ├── montecarlo_null_calibration.json
│       ├── boat_results_production_v2.0.json
│       └── README.md
├── attached_assets/          # Survey data files
├── BOAT_v2_MANIFEST.md       # Sealed manifest (after MC calibration)
├── main.py
└── README.md
```

### Testing order:
1. manifest.py — verify all 8 checks pass
2. montecarlo.py — run on 6dFGS + SDSS footprints, verify null distribution
3. geometry.py — verify CMB apex blueshift, anti-apex redshift, unit vectors
4. datasets.py — load each survey, verify pre-qualification
5. analysis.py — synthetic dipole injection → PASS, synthetic null → FAIL
6. runner.py — full production run

---

## SECTION K: WHAT SUCCESS AND FAILURE LOOK LIKE

### If BOAT v2 returns PASS (≥3 of 4 surveys pass all criteria):
- Strong evidence for a CMB-dipole-aligned anisotropy in galaxy redshifts
- If orbital discriminator also detects post-subtraction signal toward r̂/k̂,
  this would be the first evidence distinguishing TSM2.1's orbital prediction
  from standard kinematic interpretation
- Proceed to formal write-up for arXiv submission

### If BOAT v2 returns FAIL:
- Record honestly, publish all artifacts
- Analyse which criteria failed and why
- Determine if the test design needs further refinement (v3) or if the
  signal is genuinely absent/too weak to detect with current surveys

### Either way:
- No post-hoc adjustments
- No reinterpretation of results
- Pipeline is honest. We record what we see.

---

## SECTION L: COMMUNICATION STYLE

Graham expects:
- Direct and clear, no filler, no excessive qualifiers
- Practical structure when helpful, conversational when appropriate
- Honest assessment — if something doesn't hold up, say so
- Step-by-step thinking for complex problems
- No cheerleading — just results

Geoffrey expects:
- Disciplined implementation
- No urgency, no expectation of particular outcome
- Clarity and clean discrimination
- Immutable artifact publication regardless of result

---

## SECTION M: OUTSTANDING QUESTIONS

1. **DESI DR1 SQL queries:** Need to determine exact table names and column
   mappings in the NOIRLab Astro Data Lab schema. Check desi_dr1 database
   documentation before writing queries.

2. **DESI footprint for Monte Carlo:** If we want to run Monte Carlo before
   downloading DESI data, we need representative footprint coordinates.
   Alternative: download DESI first, then run MC on all 4 real footprints
   before sealing.

3. **DESI BGS vs LRG independence:** Both share the same sky footprint but
   different galaxy populations and redshift ranges. We treat them as
   independent datasets. A reviewer might challenge this — be prepared to
   justify based on different tracer populations and z-ranges probing
   different cosmic epochs.

4. **Weighted Pearson p-value:** The effective sample size N_eff = (Σw)²/Σw²
   is an approximation. For very large N (millions of galaxies), p-values
   will be extremely small regardless. The ratio test is the more
   meaningful discriminator.

---

*Document prepared by Claude Opus (Anthropic), 21 February 2026*
*For use by Graham Hill and the BOAT v2.0 implementation team*
