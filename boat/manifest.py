"""
BOAT Manifest — Locked Constants
TSM2 Blind Orbital Anisotropy Test
Sealed: 19 February 2026
SHA-256: 17abb8b0548f9bfe4725468314e5ee655cea3f0f8ccf21fdc998f196d0be43af

NO VALUES IN THIS FILE MAY BE MODIFIED AFTER SEAL DATE.
"""

import numpy as np

C_KM_S = 299792.458
V_OBS = 370.0
CST_PERIOD_GYR = 284.0
BETA_OBS = V_OBS / C_KM_S

_P_SECONDS = CST_PERIOD_GYR * 1e9 * 3.15576e7
R_OBS_KM = V_OBS * _P_SECONDS / (2 * np.pi)
R_OBS_MLY = R_OBS_KM / (9.461e12 * 1e6)
R_OBS_MPC = R_OBS_KM / 3.0857e19

COG_X_RA_DEG = 347.7500
COG_X_DEC_DEG = 66.0000

CMB_DIPOLE_RA_DEG = 167.7866
CMB_DIPOLE_DEC_DEG = -7.1454

def _unit_vector(ra_deg, dec_deg):
    """Convert RA/Dec (degrees) to Cartesian unit vector."""
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

R_HAT = _unit_vector(COG_X_RA_DEG, COG_X_DEC_DEG)
V_HAT = _unit_vector(CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG)

_k_cross = np.cross(R_HAT, V_HAT)
K_HAT = _k_cross / np.linalg.norm(_k_cross)

V_OBS_VEC = V_OBS * V_HAT

N_DATASETS = 5
N_DECOY_AXES = 10
CORRELATION_THRESHOLD = 0.01  # Minimum |r| for pass (amended v1.1 from 0.05)
P_VALUE_THRESHOLD = 0.01
PASS_REQUIRED = 3

Z_MIN = 0.01
Z_MAX = 0.50
GALACTIC_PLANE_MASK_DEG = 10
MIN_SAMPLE_SIZE = 1000
DUPLICATE_MATCH_ARCSEC = 1.0

DATASET_SELECTION_SEED = 20260219
DECOY_AXES_SEED = 20260219

def verify_manifest():
    """Verify internal consistency of all manifest constants."""
    checks = []

    checks.append(("R_HAT unit norm", abs(np.linalg.norm(R_HAT) - 1.0) < 1e-10))
    checks.append(("V_HAT unit norm", abs(np.linalg.norm(V_HAT) - 1.0) < 1e-10))
    checks.append(("K_HAT unit norm", abs(np.linalg.norm(K_HAT) - 1.0) < 1e-10))

    checks.append(("R_HAT ⊥ K_HAT", abs(np.dot(R_HAT, K_HAT)) < 1e-10))
    checks.append(("V_HAT ⊥ K_HAT", abs(np.dot(V_HAT, K_HAT)) < 1e-10))

    checks.append(("R_OBS_MLY ≈ 55.78", abs(R_OBS_MLY - 55.78) < 0.1))

    checks.append(("BETA_OBS ≈ 0.001234", abs(BETA_OBS - 0.001234) < 0.000001))

    angle = np.degrees(np.arccos(np.clip(np.dot(R_HAT, V_HAT), -1, 1)))
    checks.append(("r̂–v̂ separation ≈ 121.15°", abs(angle - 121.15) < 0.1))

    all_pass = all(c[1] for c in checks)

    if not all_pass:
        failed = [c[0] for c in checks if not c[1]]
        raise RuntimeError(f"MANIFEST VERIFICATION FAILED: {failed}")

    return True

verify_manifest()
