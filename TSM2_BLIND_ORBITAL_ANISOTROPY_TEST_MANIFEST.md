# TSM2.1 BLIND ORBITAL ANISOTROPY TEST — PRE-REGISTERED MANIFEST

**Document ID:** TSM2-BOAT-MANIFEST-v1.1  
**Date:** 19 February 2026  
**Status:** SEALED — No modifications permitted after SHA-256 hash is published  
**Authority:** TSM2 Institute for Cosmology Ltd  
**Repository:** github.com/TSM2Institute  

---

## 1. PURPOSE

This manifest pre-registers all constants, geometry, dataset rules, statistical
thresholds, and output specifications for the Blind Orbital Anisotropy Test (BOAT).

The test asks: **Does the CoG X orbital geometry produce a statistically significant
anisotropy in redshift residuals compared to random sky directions?**

This is a falsifiable, binary test. PASS or FAIL. No interpretation layer.

---

## 2. LOCKED CONSTANTS

### 2.1 Primary Constants (Observationally Anchored)

| Parameter | Value | Source | Status |
|-----------|-------|--------|--------|
| v_obs | 370 km/s | CMB dipole (Planck 2018) | Observational — non-negotiable |
| CST_PERIOD_GYR | 284.0 ± 2 Gyr | TSM2.1 (92.5 Gyr base × 3.07 UTS stretch) | Locked in repo config.py |
| c | 299,792.458 km/s | SI definition | Exact |

### 2.2 Derived Constants

| Parameter | Value | Derivation |
|-----------|-------|------------|
| R_obs | 55.78 Mly (17.10 Mpc) | v_obs × P / (2π) |
| β_obs | 0.001234 | v_obs / c |
| ω | 2π / P | Angular frequency |

### 2.3 Phase 1 Pipeline Constants (Reference Only — Not Used in BOAT Prediction)

| Parameter | Value | Notes |
|-----------|-------|-------|
| K_TSM | 5.1 × 10⁻²³ cm² | Scattering coefficient |
| N_COSMIC_BASELINE_HIGHZ | 2.5 × 10²⁰ cm⁻² | HI column baseline |
| COSMIC_EXPONENT | 2.3 | Power-law index |
| B_FIELD | 10⁻⁶ Gauss | Intergalactic magnetic field |

These constants are documented for traceability. They are NOT inputs to the
BOAT z_model prediction, which uses only orbital geometry.

---

## 3. COORDINATE SYSTEM AND GEOMETRY

### 3.1 Coordinate Frame

All coordinates expressed in **ICRS (International Celestial Reference System)**,
epoch J2000.0. Unit vectors computed from (RA, Dec) via standard spherical-to-Cartesian
conversion:

```
n̂ = [ cos(Dec) × cos(RA),  cos(Dec) × sin(RA),  sin(Dec) ]
```

where RA and Dec are in radians.

### 3.2 Primary Directional Vectors

**CoG X radial direction (r̂):**
- RA = 347.7500° (23h 11m 00s), Dec = +66.0000°
- Galactic: l = 112.9532°, b = +5.0934°
- Unit vector: [0.397476, −0.086300, 0.913545]

**CMB dipole velocity direction (v̂):**
- Galactic: l = 264.00°, b = +48.00° (Planck 2018)
- ICRS: RA = 167.7866°, Dec = −7.1454°
- Unit vector: [−0.969776, 0.209910, −0.124388]

**Angular separation r̂ vs v̂:** 121.15°

Note: The 31.15° departure from the 90° expected for a purely circular orbit is
acknowledged. The test does not assume circularity — it anchors velocity to the
observational CMB dipole and tests whether the CoG X centreline produces a
reproducible anisotropy signature.

### 3.3 Derived Orbital Plane Normal (k̂)

Computed as the angular-momentum direction:

```
k̂ = (r̂ × v̂) / |r̂ × v̂|
```

- ICRS: RA = 257.7888°, Dec = −0.0173°
- Galactic: l = 20.8742°, b = +22.1469°
- Unit vector: [−0.211516, −0.977374, −0.000301]

**Verification (all confirmed computationally):**
- r̂ · k̂ = 0.000000 (orthogonal)
- v̂ · k̂ = 0.000000 (orthogonal)
- |k̂| = 1.000000 (normalised)

### 3.4 Observer Velocity Vector

```
v⃗_obs = 370.0 × v̂ = [−358.8171, +77.6668, −46.0234] km/s
```

### 3.5 Phase Angle (φ)

With v̂ anchored directly from the CMB dipole, the phase angle φ is implicit
in the geometry and is not an independent parameter. It is omitted from the
manifest to avoid redundancy.

---

## 4. PREDICTION METHOD

### 4.1 z_model Computation (Redshift-Independent)

For each galaxy with sky direction n̂ (from RA, Dec):

**Step 1.** Convert RA/Dec to unit vector:
```
n̂_gal = [ cos(Dec) × cos(RA),  cos(Dec) × sin(RA),  sin(Dec) ]
```

**Step 2.** Compute projected line-of-sight velocity:
```
v_LOS = −v⃗_obs · n̂_gal
```
(Negative sign: motion toward source = blueshift)

**Step 3.** Compute projected β:
```
β_proj = v_LOS / c
```

**Step 4.** Apply SR Doppler:
```
z_model = √((1 + β_proj) / (1 − β_proj)) − 1
```

**No observed redshift is used at any step in computing z_model.**

### 4.2 Expected Signal Characteristics

| Direction | β_proj | z_model |
|-----------|--------|---------|
| Toward CMB apex (max blueshift) | −0.001234 | −0.001233 |
| Toward CMB anti-apex (max redshift) | +0.001234 | +0.01235 |
| Toward CoG X | +0.000638 | +0.000638 |
| Dipole amplitude | — | ~0.00124 |

The signal is a dipole pattern aligned with the CMB dipole direction, with
potential higher-order (quadrupole) structure from orbital curvature.

---

## 5. DATASET INCLUSION RULES

### 5.1 Survey Selection

From publicly available wide-field spectroscopic surveys, randomly select
**K = 5 datasets** using a recorded random seed. Candidate surveys include
(but are not limited to):

- SDSS DR18 (Sloan Digital Sky Survey)
- 2dFGRS (2-degree Field Galaxy Redshift Survey)
- 6dFGS (6-degree Field Galaxy Survey)
- GAMA (Galaxy and Mass Assembly)
- WiggleZ Dark Energy Survey

Selection criteria: wide sky coverage (>1000 deg²), spectroscopic redshifts,
publicly accessible catalogues.

### 5.2 Quality Cuts

For each dataset, apply the following filters before analysis:

| Filter | Criterion |
|--------|-----------|
| Redshift quality | Only secure/reliable quality flags (survey-specific) |
| Redshift range | 0.01 ≤ z ≤ 0.5 (low-z to maximise orbital dipole signal-to-noise) |
| Valid coordinates | Finite RA, Dec within valid ranges |
| Duplicate removal | Unique objects only (nearest-neighbour matching at 1 arcsec) |
| Galactic plane mask | Exclude |b| < 10° (Zone of Avoidance contamination) |
| Minimum sample size | N ≥ 1000 objects per dataset after cuts |

### 5.3 Redshift Range Justification

The orbital dipole signal (~0.001 in z) is most detectable relative to observed
redshifts at low-z (z < 0.5) where it represents a larger fractional contribution.
At high-z, the dominant distance-dependent redshift (z >> 1) overwhelms the
orbital dipole signal.

### 5.4 Dataset Versioning

Record for each dataset: catalogue name, data release version, download date,
total objects before cuts, total objects after cuts. Compute SHA-256 hash of
each filtered dataset CSV and include in the output artifact.

---

## 6. STATISTICAL ANALYSIS SPECIFICATION

### 6.1 Residual Computation

For each galaxy i:
```
Δz_i = z_obs,i − z̄_obs − z_model,i
```

where z̄_obs is the sample mean observed redshift. The subtraction of z̄_obs
removes the isotropic redshift component (dominated by Hubble flow / distance),
isolating the directional dipole signal.

### 6.2 Regression Model

Fit the following model to residuals:

**Dipole only (primary):**
```
Δz = A × cos(θ) + ε
```
where θ is the angular separation between galaxy direction n̂ and velocity
direction v̂, and A is the dipole amplitude.

**Dipole + Quadrupole (pre-registered extension):**
```
Δz = A × cos(θ) + B × (3cos²(θ) − 1)/2 + ε
```

Report which term(s) are significant at p < 0.01.

### 6.3 Decoy Control Axes

For each dataset:
- Generate **N_decoy = 10** random sky axes (uniform on sphere, recorded seed)
- Run identical dipole + quadrupole analysis on each decoy axis
- Record r, p-value, RMS residual for each axis
- The true CoG X / CMB velocity axis is blinded among the 11 total axes

Do not reveal which axis is real until all 11 results are computed and recorded.

### 6.4 PASS/FAIL Thresholds

**Per-dataset PASS requires ALL of:**

| Criterion | Threshold |
|-----------|-----------|
| Correlation coefficient | \|r\| ≥ 0.01 |
| p-value | p < 0.01 |
| True axis rank | Highest \|r\| among all 11 axes |

Note: The correlation threshold is set at |r| ≥ 0.01 (revised from 0.05 in v1.0 after synthetic validation revealed the expected signal-to-background ratio of ~0.009 makes |r| ≥ 0.05 physically unreachable). The p-value and rank criteria provide the primary discrimination power. Amendment approved by Geoffrey Thwaites, Claude (Anthropic), and Grok (xAI) on 19 Feb 2026.

**Overall PASS requires:**
- ≥ 3 of 5 datasets pass all three per-dataset criteria

**Otherwise: FAIL.**

No interpretation. Only numbers.

---

## 7. OUTPUT ARTIFACTS

### 7.1 Required Outputs (Per Dataset)

| Output | Description |
|--------|-------------|
| Dataset name/version | Catalogue identifier and release |
| N_objects | Sample size after cuts |
| r_true | Correlation coefficient (true axis) |
| p_true | p-value (true axis) |
| RMS_true | RMS residual (true axis) |
| r_decoy_1 ... r_decoy_10 | Correlation coefficients (all decoy axes) |
| Rank_true | Rank of true axis among 11 (1 = highest \|r\|) |
| Significant_terms | Which regression terms significant at p < 0.01 |
| PASS / FAIL | Binary verdict per dataset |

### 7.2 Required Outputs (Overall)

| Output | Description |
|--------|-------------|
| N_pass | Count of datasets that passed |
| N_fail | Count of datasets that failed |
| Overall verdict | PASS if N_pass ≥ 3, else FAIL |

### 7.3 Immutable Artifact Package

The following files constitute the sealed test record:

1. **This manifest** (SHA-256 hashed at seal time)
2. **Code commit hash** (Git SHA of the test implementation)
3. **Dataset hashes** (SHA-256 of each filtered dataset CSV)
4. **Results JSON** (machine-readable per-dataset and overall results)
5. **Final report PDF** (human-readable summary with all tables)
6. **Random seeds file** (all seeds used for dataset selection and decoy generation)

All artifacts committed to github.com/TSM2Institute repository in a single
timestamped commit. No subsequent modifications to this commit.

---

## 8. NON-NEGOTIABLE RULES

1. No parameter tuning after manifest hash is published
2. No axis adjustment after manifest hash is published
3. No tolerance changes after manifest hash is published
4. No filtering changes after results are observed
5. No dataset substitution after random selection is performed
6. No re-running with different seeds
7. No "exploratory" pre-analysis on the selected datasets
8. Results are published regardless of outcome (PASS or FAIL)

**This is a strict blind test. Science wins clean.**

---

## 9. APPROVAL AND SEAL

**Theory:** Geoffrey E. Thwaites — TSM2.1 Author  
**Verification Systems:** Claude (Anthropic), Grok (xAI), ChatGPT (OpenAI)  
**Implementation:** Graham Hill — TSM2 Institute  

**Manifest SHA-256:** [SEE COMMIT MESSAGE]  
**Seal Date:** 19 February 2026  

---

*"The simplest test is the most decisive." — TSM2 Institute*
