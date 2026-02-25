# BOAT v2.0 — Blind Orbital Anisotropy Test

## Pre-Registered Manifest

**Document ID:** BOAT-MANIFEST-v2.0  
**Date:** [TO BE SEALED AFTER MONTE CARLO CALIBRATION]  
**Authors:** Geoffrey E. Thwaites, Graham Hill, Claude (Anthropic), Grok (xAI)  
**Institution:** TSM2 Institute for Cosmology Ltd  
**Repository:** github.com/TSM2Institute/blind-orbital-anisotropy-test  
**Predecessor:** BOAT v1.1 (sealed 19 February 2026, overall verdict FAIL)

---

## 1. Objective

Test whether galaxy redshift residuals exhibit a statistically significant
anisotropy aligned with the CMB dipole velocity direction, consistent with
the TSM2.1 orbital motion prediction, and whether post-dipole residuals
show additional structure aligned with the predicted Centre of Gravity
direction (CoG X).

This test is executed blind: all constants, thresholds, axes, and datasets
are locked before any real data touches the analysis pipeline. The result
is published regardless of outcome.

---

## 2. Changes from v1.1

| Feature | v1.1 | v2.0 |
|---------|------|------|
| PASS criterion | Rank must be #1 of 11 | Ratio ≥ R_NULL_99 AND top 3 of 51 |
| Decoy axes | 10 (seeded) | 50 (seeded) |
| Survey selection | Any with N ≥ 1000 | Pre-qualified: sky ≥ 15%, median z ≥ 0.04 |
| Footprint correction | None | Inverse local density weighting (R = 5°) |
| z-bin invariance | Not tested | 4-bin consistency required |
| Orbital discriminator | CMB dipole only | CMB dipole + post-subtraction CoG X / k̂ |
| Overall PASS rule | 3 of 5 | ⌈N/2⌉ + 1 of N (expected: 3 of 4) |
| Ratio threshold | Heuristic | Monte Carlo calibrated (global FPR < 0.01) |

---

## 3. Locked Physical Constants

All constants are inherited from v1.1. No changes.

### 3.1 Observational Anchors

| Constant | Value | Source |
|----------|-------|--------|
| v_obs | 370 km/s | CMB dipole (Planck 2018) |
| CMB dipole apex | RA 167.7866°, Dec −7.1454° (J2000) | Planck 2018 |
| CMB dipole anti-apex | RA 347.7866°, Dec +7.1454° (J2000) | Derived |

### 3.2 TSM2.1 Orbital Parameters

| Constant | Value | Derivation |
|----------|-------|------------|
| CST_PERIOD_GYR | 284.0 ± 2 | Repository locked |
| R_obs | 55.78 Mly | v_obs × P / (2π) |
| CoG X direction | RA 347.75°, Dec +66.00° (J2000) | Geoffrey Thwaites (2018) |

### 3.3 Derived Vectors (ICRS J2000.0)

| Vector | Cartesian [x, y, z] | Definition |
|--------|---------------------|------------|
| r̂ (CoG X) | [+0.397476, −0.086300, +0.913545] | Radial toward centre |
| v̂ (CMB dipole) | [−0.969776, +0.209910, −0.124388] | Velocity direction |
| k̂ (orbital normal) | [−0.211516, −0.977374, −0.000301] | r̂ × v̂ / |r̂ × v̂| |
| v⃗_obs (km/s) | [−358.8171, +77.6668, −46.0234] | v_obs × v̂ |

### 3.4 Verification Checks (run on import)

1. |r̂| = 1.000000 ± 10⁻⁶
2. |v̂| = 1.000000 ± 10⁻⁶
3. |k̂| = 1.000000 ± 10⁻⁶
4. r̂ · k̂ = 0.000000 ± 10⁻⁶
5. v̂ · k̂ = 0.000000 ± 10⁻⁶
6. R_obs ≈ 55.78 Mly ± 0.01
7. β_obs ≈ 0.001234 ± 10⁻⁶
8. Angular separation r̂–v̂ ≈ 121.15° ± 0.01°

---

## 4. Redshift Model (Unchanged from v1.1)

The predicted redshift contribution from orbital motion:

1. Convert galaxy (RA, Dec) → unit vector n̂
2. Compute line-of-sight velocity: v_LOS = −v⃗_obs · n̂
3. Compute β = v_LOS / c
4. Apply SR Doppler: z_model = √((1+β)/(1−β)) − 1

**No observed redshift is used in computing z_model.**

---

## 5. Survey Pre-Qualification

### 5.1 Qualification Criteria

A survey qualifies for BOAT v2.0 if and only if:

| Criterion | Threshold |
|-----------|-----------|
| Sky coverage fraction | ≥ 15% of full sphere |
| Median spectroscopic redshift | ≥ 0.04 |
| Minimum sample after quality cuts | ≥ 10,000 galaxies |
| Spectroscopic (not photometric) redshifts | Required |
| Public data access | Required |

Sky fraction is computed as: number of occupied 1° × 1° (RA, Dec) bins
divided by 64,800 (total bins on sphere).

### 5.2 Pre-Qualified Surveys

| # | Survey | Expected Sky % | Expected Median z | Hemisphere | Access |
|---|--------|---------------|-------------------|------------|--------|
| 1 | 6dFGS DR3 | ~34% | ~0.05 | Southern | VizieR/ROE |
| 2 | SDSS DR18 | ~21% | ~0.20 | Northern | SciServer CasJobs |
| 3 | DESI DR1 BGS | ~18% | ~0.20 | Both (northern-heavy) | NERSC |
| 4 | DESI DR1 LRG | ~18% | ~0.70 | Both (northern-heavy) | NERSC |

DESI BGS and LRG are treated as independent datasets because they probe
different redshift regimes (z < 0.6 vs 0.4 < z < 1.0) with different
galaxy populations, despite sharing the same sky footprint.

If any survey fails pre-qualification after data loading (e.g., sky fraction
falls below 15% after quality cuts), it is excluded and the overall PASS
rule adjusts to ⌈(N-1)/2⌉ + 1.

### 5.3 Quality Cuts (Applied to All Surveys)

| Filter | Criterion |
|--------|-----------|
| Redshift range | 0.01 ≤ z_obs ≤ 0.5 (BGS, 6dFGS, SDSS) or 0.4 ≤ z_obs ≤ 1.0 (LRG) |
| Galactic plane mask | \|b\| ≥ 10° |
| Duplicate removal | Nearest-neighbour < 1 arcsec (astropy match_coordinates_sky) |
| Minimum sample | N ≥ 10,000 after all cuts |
| Survey-specific quality flags | Per-survey reliable redshift criteria |

---

## 6. Footprint Correction

### 6.1 Method

Inverse local density weighting compensates for non-uniform sky coverage.

For each galaxy i:
1. Count N_local = number of galaxies within angular radius R_weight
2. Compute weight w_i = 1 / N_local
3. Normalise weights: w_i → w_i / mean(w_i)

Weighted Pearson correlation and weighted linear regression use these
weights throughout all statistical analyses.

### 6.2 Weighting Radius

**R_weight = 5° (fixed)**

This is the standard two-halo transition scale in galaxy angular clustering.
It is locked without calibration on any real or mock survey data.
No optimisation of this parameter is permitted.

---

## 7. Statistical Analysis

### 7.1 Residual Computation

Δz_i = z_obs,i − mean(z_obs)

The sample mean is subtracted to remove the isotropic Hubble-flow
component. z_model is computed for reference but NOT subtracted from
the data.

### 7.2 Dipole Regression (Primary Test)

For each test axis a⃗:

1. Compute cos(θ_i) = n̂_i · a⃗ for all galaxies
2. Fit weighted linear model: Δz = A × cos(θ) + intercept
3. Compute weighted Pearson correlation r_weighted
4. Compute p-value from weighted correlation

### 7.3 Dipole + Quadrupole Regression

Δz = A × cos(θ) + B × (3cos²(θ) − 1)/2 + intercept

Report both dipole (A) and quadrupole (B) amplitudes with p-values.

### 7.4 Decoy Axes

- N_decoy = 50 (increased from 10 in v1.1)
- Generated uniform on sphere from seeded RNG
- Seed: DECOY_AXES_SEED = 20260301
- Total axes tested per dataset: 51 (1 true + 50 decoys)

### 7.5 PASS Criteria (Per Dataset)

A dataset PASSES if ALL of the following hold:

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Ratio test | true \|r\| / mean(decoy \|r\|) ≥ R_NULL_99 | Exceeds 99th percentile of null distribution |
| Rank test | True axis in top 3 of 51 axes | Safety gate against outlier decoys |
| Significance | p < 0.01 on true axis | Statistical significance |
| z-bin invariance | max/min \|r\| ratio ≤ 2.0 across qualifying bins | Signal is distance-independent |

**R_NULL_99** is the ratio threshold corresponding to a global
false-positive rate < 0.01 under the null hypothesis. It is determined
by Monte Carlo simulation (Section 7.8) and locked before any real data
is analysed.

### 7.6 Redshift-Binned Invariance Test

Split each dataset into fixed bins:

**Standard bins (6dFGS, SDSS, DESI BGS):**

| Bin | Range |
|-----|-------|
| 1 | 0.01 – 0.05 |
| 2 | 0.05 – 0.10 |
| 3 | 0.10 – 0.20 |
| 4 | 0.20 – 0.50 |

**LRG bins (DESI LRG):**

| Bin | Range |
|-----|-------|
| 1 | 0.40 – 0.55 |
| 2 | 0.55 – 0.70 |
| 3 | 0.70 – 0.85 |
| 4 | 0.85 – 1.00 |

Rules:
- Bins with < 1,000 galaxies are excluded from invariance test
- At least 2 qualifying bins required; otherwise invariance test is waived
- True axis |r| must have the same sign across all qualifying bins
- max(|r|) / min(|r|) ≤ 2.0 across qualifying bins

### 7.7 Weighted Statistics Implementation

All correlations and regressions use the footprint weights from Section 6.

Weighted Pearson correlation:

r_w = Σ w_i (x_i − x̄_w)(y_i − ȳ_w) / √[Σ w_i (x_i − x̄_w)² × Σ w_i (y_i − ȳ_w)²]

where x̄_w and ȳ_w are weighted means.

Effective sample size for p-value computation:

N_eff = (Σ w_i)² / Σ w_i²

P-value computed from t-distribution with N_eff − 2 degrees of freedom.

### 7.8 Monte Carlo Null Calibration

**Purpose:** Determine R_NULL_99, the ratio threshold that yields a global
false-positive rate < 0.01 under the null hypothesis (no dipole signal).

**Method:**

For each of the 4 pre-qualified surveys:

1. Generate synthetic isotropic mock catalogues matching the survey's:
   - Sky footprint (same RA/Dec distribution as real data)
   - Galaxy count (same N after quality cuts)
   - Redshift distribution (same z histogram as real data)
   - No injected dipole signal (z_obs drawn from the real z distribution
     but randomly reassigned to positions)
2. Run the complete v2 analysis protocol on each mock:
   - Apply footprint weighting (R = 5°)
   - Compute true axis and 50 decoy axes
   - Record ratio = true |r| / mean(decoy |r|)
3. Repeat for 10,000 independent realisations per survey
4. Combine all 40,000 ratio values (4 surveys × 10,000 mocks)
5. R_NULL_99 = 99th percentile of the combined null ratio distribution

**Locked value:** R_NULL_99 = [TO BE FILLED BY MONTE CARLO — locked before
any real v2 analysis]

**Documentation:** The full distribution of null ratios, per-survey
breakdowns, and the 95th/99th/99.5th percentile values are recorded
in the Monte Carlo results file (SHA-256 hashed and committed).

---

## 8. Orbital Discriminator (New in v2.0)

### 8.1 Purpose and Physical Basis

Standard cosmology predicts a kinematic dipole aligned with the CMB
dipole direction (v̂), arising from our peculiar motion relative to the
CMB rest frame. This is well-established and uncontroversial.

TSM2.1 predicts this same dipole but attributes it to orbital motion
around a Centre of Gravity (CoG X). The orbital interpretation carries
a specific, testable consequence: the velocity vector v̂ and the radial
direction r̂ (toward CoG X) are separated by 121.15° — not the 90°
expected for a circular orbit. This offset implies a non-circular
(elliptical or spiral) orbital trajectory with a phase-dependent
relationship between the radial and tangential velocity components.

After subtracting the dominant v̂ dipole from the data, any residual
anisotropy aligned with r̂ (radial toward CoG X) or k̂ (orbital normal,
perpendicular to the orbital plane) would represent structure that
standard kinematic models do not predict. Specifically:

- A residual toward r̂ indicates a radial velocity component consistent
  with non-circular orbital motion (approach toward or recession from
  CoG X beyond what the v̂ projection accounts for).
- A residual toward k̂ indicates structure perpendicular to the orbital
  plane, which could arise from orbital precession or out-of-plane
  perturbations.

This discriminator is hypothesis-driven by the 121° r̂–v̂ separation
and the non-circular orbital geometry. It is not exploratory.

### 8.2 Method

**Step 1: Primary dipole fit**

Fit Δz = A_v × cos(θ_v) + intercept, where θ_v is the angle from v̂.

**Step 2: Subtract fitted dipole**

Δz_residual = Δz − A_v × cos(θ_v) − intercept

**Step 3: Test residual against r̂ and k̂**

For each axis (r̂ and k̂):
1. Compute cos(θ) from the axis
2. Fit weighted model: Δz_residual = A × cos(θ) + intercept
3. Compute weighted Pearson r and p-value

**Step 4: Compare against decoys**

Run the same post-subtraction test on 50 decoy axes.
Compute ratio = |r_axis| / mean(|r_decoy|) for both r̂ and k̂.

### 8.3 Orbital Discriminator PASS Criteria

The orbital discriminator is a SECONDARY test. It does not affect the
primary PASS/FAIL verdict. It is reported separately.

| Criterion | Threshold |
|-----------|-----------|
| Post-subtraction p-value | < 0.05 on r̂ or k̂ |
| Post-subtraction ratio | \|r\| / mean(decoy \|r\|) ≥ 1.5 |

A **DETECTION** is reported if both criteria are met for r̂ or k̂.
A **NON-DETECTION** is reported otherwise.

---

## 9. Overall Verdict

### 9.1 Primary Verdict

For N qualifying surveys (expected N = 4):

**PASS** if ≥ ⌈N/2⌉ + 1 datasets pass ALL primary criteria (Section 7.5).

| N | Required passes |
|---|----------------|
| 3 | 3 (unanimous) |
| 4 | 3 |
| 5 | 4 |

**FAIL** otherwise.

### 9.2 Orbital Discriminator Summary

Reported separately as:
- "DETECTED in K of N surveys" (if any show post-subtraction signal)
- "NOT DETECTED in any survey"

---

## 10. Immutable Artifacts

All artifacts are committed to the TSM2Institute GitHub repository
in a single timestamped commit after execution is complete.

| Artifact | Format | Hash |
|----------|--------|------|
| This manifest | .md | SHA-256 (sealed before execution) |
| Monte Carlo null calibration results | .json | SHA-256 |
| All source code | .py | Git commit SHA |
| Decoy axes file | .json | SHA-256 |
| Per-dataset filtered data hashes | SHA-256 | In results JSON |
| Per-dataset results | .json | SHA-256 |
| Production results (complete) | .json | SHA-256 |
| Final report | .md | SHA-256 |

---

## 11. Execution Protocol

1. Build Monte Carlo calibration module
2. Run null simulation (10,000 × 4 surveys) → determine R_NULL_99
3. Fill R_NULL_99 in this manifest (Section 7.8)
4. Seal this manifest (compute SHA-256, commit to GitHub)
5. Download and load all 4 survey catalogues
6. Verify pre-qualification (sky ≥ 15%, median z ≥ 0.04)
7. Execute blind test on each qualifying survey
8. Run orbital discriminator on each survey
9. Compile results
10. Publish all artifacts regardless of outcome

No modifications to any parameter, threshold, or criterion after Step 4.

---

## 12. Signatories

**Theory:** Geoffrey E. Thwaites — TSM2.1 Author  
**Implementation & Coordination:** Graham Hill — TSM2 Institute  
**Verification:** Claude (Anthropic), Grok (xAI)

---

## 13. SHA-256 Hash

**Manifest hash (pre-seal):** [COMPUTED AT SEAL TIME, AFTER R_NULL_99 IS FILLED]

---

*"The simplest test is the most decisive."*

*Pipeline is honest. We record what we see.*
