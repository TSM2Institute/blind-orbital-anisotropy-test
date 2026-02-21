# BOAT v2.0 — Build Prompt #2: Data File Rename + manifest.py

**Context:** The clone has been cleaned (prompt #1 complete). geometry.py and hashing.py are kept. geometry.py cannot import yet because it depends on manifest.py. This prompt fixes that.

---

## TASK 1: RENAME DATA FILES

The data files in `attached_assets/` have Replit upload timestamp suffixes. Rename them to clean, predictable names that all downstream modules will reference.

Run this once:

```python
import os
import shutil

asset_dir = "attached_assets"
renames = {
    "6dFGSzDR3.txt_1771500188269.gz": "6dFGS_DR3.txt.gz",
    "MyTable_TSM2_0_1771642595705.csv": "SDSS_DR18.csv"
}

for old, new in renames.items():
    old_path = os.path.join(asset_dir, old)
    new_path = os.path.join(asset_dir, new)
    if os.path.exists(old_path) and not os.path.exists(new_path):
        shutil.move(old_path, new_path)
        print(f"Renamed: {old} → {new}")
    elif os.path.exists(new_path):
        print(f"Already clean: {new}")
    else:
        print(f"WARNING: {old} not found!")
```

After running, confirm both files exist with their new names.

**Do not rename any other files.** The unused survey files (2dFGRS, GAMA, 2MRS) keep their timestamped names.

---

## TASK 2: BUILD boat/manifest.py

Create the file `boat/manifest.py` with the exact content below. Do not modify any constant values — these are locked and verified across multiple independent systems.

```python
"""
BOAT v2.0 — Locked Manifest Constants

All physical constants, derived vectors, and test parameters for the
Blind Orbital Anisotropy Test v2.0.

NOTHING IN THIS FILE MAY BE CHANGED AFTER THE MANIFEST IS SEALED.
The only exception is R_NULL_99, which is filled by Monte Carlo
calibration BEFORE sealing.

Verification checks run automatically on import.
"""

import numpy as np
import math

# ============================================================
# OBSERVATIONAL ANCHORS
# ============================================================

V_OBS = 370.0                # km/s — CMB dipole speed (Planck 2018)
C_KM_S = 299792.458          # km/s — speed of light

CMB_DIPOLE_RA_DEG = 167.7866     # degrees (J2000)
CMB_DIPOLE_DEC_DEG = -7.1454     # degrees (J2000)

# Anti-apex (derived)
CMB_ANTI_APEX_RA_DEG = 347.7866  # degrees (J2000)
CMB_ANTI_APEX_DEC_DEG = 7.1454   # degrees (J2000)

# ============================================================
# TSM2.1 ORBITAL PARAMETERS
# ============================================================

CST_PERIOD_GYR = 284.0       # Gyr — orbital period
COG_X_RA_DEG = 347.75        # degrees (J2000) — Centre of Gravity direction
COG_X_DEC_DEG = 66.00        # degrees (J2000)

# ============================================================
# DERIVED CONSTANTS
# ============================================================

R_OBS_MLY = 55.78             # Mly — orbital radius (V_OBS × P / 2π)
BETA_OBS = V_OBS / C_KM_S    # ≈ 0.001234

# ============================================================
# DERIVED UNIT VECTORS (ICRS J2000.0)
# ============================================================

def _ra_dec_to_unit(ra_deg, dec_deg):
    """Convert RA, Dec (degrees) to Cartesian unit vector."""
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

# r̂ — radial direction toward Centre of Gravity (CoG X)
R_HAT = _ra_dec_to_unit(COG_X_RA_DEG, COG_X_DEC_DEG)
# Expected: [+0.397476, -0.086300, +0.913545]

# v̂ — velocity direction (CMB dipole apex)
V_HAT = _ra_dec_to_unit(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)
# Expected: [-0.969776, +0.209910, -0.124388]

# k̂ — orbital normal: r̂ × v̂ / |r̂ × v̂|
_k_raw = np.cross(R_HAT, V_HAT)
K_HAT = _k_raw / np.linalg.norm(_k_raw)
# Expected: [-0.211516, -0.977374, -0.000301]

# Velocity vector (km/s)
V_OBS_VEC = V_OBS * V_HAT
# Expected: [-358.8171, +77.6668, -46.0234]

# ============================================================
# v2.0 TEST PARAMETERS
# ============================================================

N_DECOY_AXES = 50
DECOY_AXES_SEED = 20260301
P_VALUE_THRESHOLD = 0.01
RANK_THRESHOLD = 3              # True axis must be in top 3 of 51

# Redshift cuts — standard surveys (6dFGS, SDSS, DESI BGS)
Z_MIN = 0.01
Z_MAX = 0.50

# Redshift cuts — LRG surveys (DESI LRG)
Z_MIN_LRG = 0.40
Z_MAX_LRG = 1.00

# Quality cuts
GALACTIC_PLANE_MASK_DEG = 10.0
MIN_SAMPLE_SIZE = 10000
DUPLICATE_MATCH_ARCSEC = 1.0    # arcsec, for astropy match_coordinates_sky

# Pre-qualification thresholds
MIN_SKY_FRACTION = 0.15
MIN_MEDIAN_Z = 0.04

# Footprint weighting
R_WEIGHT_DEG = 5.0              # Fixed. No calibration permitted.

# z-bin invariance
MAX_MIN_RATIO_ZBIN = 2.0
MIN_BIN_SIZE = 1000

# Standard z-bins (6dFGS, SDSS, DESI BGS)
Z_BINS_STANDARD = [
    (0.01, 0.05),
    (0.05, 0.10),
    (0.10, 0.20),
    (0.20, 0.50),
]

# LRG z-bins (DESI LRG)
Z_BINS_LRG = [
    (0.40, 0.55),
    (0.55, 0.70),
    (0.70, 0.85),
    (0.85, 1.00),
]

# Overall verdict
PASS_REQUIRED_FORMULA = "ceil(N/2) + 1"

# Orbital discriminator (secondary test)
ORBITAL_DISC_P_THRESHOLD = 0.05
ORBITAL_DISC_RATIO_THRESHOLD = 1.5

# Monte Carlo calibration
N_MC_REALISATIONS = 10000
MC_FPR_TARGET = 0.01           # Global false-positive rate target

# R_NULL_99 — Monte Carlo calibrated ratio threshold
# [TO BE FILLED AFTER MONTE CARLO SIMULATION]
# This is the ONLY value that changes before sealing.
R_NULL_99 = None

# ============================================================
# DATA FILE PATHS (clean names, relative to project root)
# ============================================================

DATA_DIR = "attached_assets"
DATA_FILES = {
    "6dFGS": "6dFGS_DR3.txt.gz",
    "SDSS": "SDSS_DR18.csv",
    # DESI files will be added after download
    # "DESI_BGS": "DESI_DR1_BGS.csv",
    # "DESI_LRG": "DESI_DR1_LRG.csv",
}

# ============================================================
# VERIFICATION CHECKS (run on import)
# ============================================================

def _verify_manifest():
    """Run all 8 verification checks. Raises AssertionError on failure."""
    tol_unit = 1e-6
    tol_angle = 0.01  # degrees
    tol_r_obs = 0.01  # Mly

    # Check 1: |r̂| = 1
    assert abs(np.linalg.norm(R_HAT) - 1.0) < tol_unit, \
        f"FAIL: |r̂| = {np.linalg.norm(R_HAT)}, expected 1.0"

    # Check 2: |v̂| = 1
    assert abs(np.linalg.norm(V_HAT) - 1.0) < tol_unit, \
        f"FAIL: |v̂| = {np.linalg.norm(V_HAT)}, expected 1.0"

    # Check 3: |k̂| = 1
    assert abs(np.linalg.norm(K_HAT) - 1.0) < tol_unit, \
        f"FAIL: |k̂| = {np.linalg.norm(K_HAT)}, expected 1.0"

    # Check 4: r̂ · k̂ = 0 (perpendicular)
    assert abs(np.dot(R_HAT, K_HAT)) < tol_unit, \
        f"FAIL: r̂ · k̂ = {np.dot(R_HAT, K_HAT)}, expected 0.0"

    # Check 5: v̂ · k̂ = 0 (perpendicular)
    assert abs(np.dot(V_HAT, K_HAT)) < tol_unit, \
        f"FAIL: v̂ · k̂ = {np.dot(V_HAT, K_HAT)}, expected 0.0"

    # Check 6: R_obs ≈ 55.78 Mly
    assert abs(R_OBS_MLY - 55.78) < tol_r_obs, \
        f"FAIL: R_obs = {R_OBS_MLY}, expected 55.78"

    # Check 7: β_obs ≈ 0.001234
    assert abs(BETA_OBS - 0.001234) < tol_unit, \
        f"FAIL: β_obs = {BETA_OBS}, expected ~0.001234"

    # Check 8: Angular separation r̂–v̂ ≈ 121.15°
    cos_angle = np.dot(R_HAT, V_HAT)
    angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
    assert abs(angle_deg - 121.15) < tol_angle, \
        f"FAIL: r̂–v̂ angle = {angle_deg:.2f}°, expected 121.15°"

    print("manifest.py: All 8 verification checks PASSED")
    print(f"  |r̂| = {np.linalg.norm(R_HAT):.6f}")
    print(f"  |v̂| = {np.linalg.norm(V_HAT):.6f}")
    print(f"  |k̂| = {np.linalg.norm(K_HAT):.6f}")
    print(f"  r̂ · k̂ = {np.dot(R_HAT, K_HAT):.6e}")
    print(f"  v̂ · k̂ = {np.dot(V_HAT, K_HAT):.6e}")
    print(f"  R_obs = {R_OBS_MLY} Mly")
    print(f"  β_obs = {BETA_OBS:.6f}")
    print(f"  r̂–v̂ angle = {angle_deg:.2f}°")

# Run verification on import
_verify_manifest()
```

**IMPORTANT:** Copy this exactly. Do not reformat, reorder, or "improve" anything. Every constant value has been cross-verified by multiple independent systems.

---

## TASK 3: VERIFY

After creating manifest.py, run these checks:

### Check 1: manifest.py imports and all 8 checks pass
```python
import boat.manifest
# Should print "All 8 verification checks PASSED" with values
```

### Check 2: geometry.py now imports successfully
```python
from boat.geometry import radec_to_unit_vector
v = radec_to_unit_vector(167.7866, -7.1454)
print(f"CMB dipole unit vector: {v}")
# Should be approximately [-0.9698, +0.2099, -0.1244]
```

### Check 3: Confirm data files renamed
```python
import os
print(os.path.exists("attached_assets/6dFGS_DR3.txt.gz"))   # True
print(os.path.exists("attached_assets/SDSS_DR18.csv"))        # True
```

---

## REPORT BACK

Once complete, report:
1. Data file renames: confirmed / any issues
2. manifest.py verification: all 8 checks PASS / any failures (with values)
3. geometry.py import: PASS / FAIL
4. Data files confirmed at clean paths: yes / no

**Do not proceed to the next module.** Wait for the next build prompt from Graham.

---

*Build prompt prepared by Claude Opus (Anthropic), 21 February 2026*
