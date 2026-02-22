# BOAT v2.0 — Build Prompt #4: analysis.py

**Context:** manifest.py verified. montecarlo.py built and validated (20,000 realisations complete, provisional R_NULL_99 = 1.5627). geometry.py and hashing.py unchanged from v1.1. We now build the production analysis engine.

---

## WHAT THIS MODULE DOES

analysis.py implements the complete BOAT v2.0 statistical analysis protocol:

1. **Footprint correction** (manifest Section 6) — inverse local density weighting
2. **Residual computation** (Section 7.1) — Δz = z_obs − mean(z_obs)
3. **Dipole regression** (Section 7.2) — weighted Pearson r + p-value for each axis
4. **Dipole + quadrupole regression** (Section 7.3) — report amplitudes
5. **Decoy axis generation** (Section 7.4) — 50 seeded decoys + 1 true
6. **PASS criteria evaluation** (Section 7.5) — ratio, rank, p-value, z-bin invariance
7. **z-bin invariance test** (Section 7.6) — 4-bin consistency check
8. **Weighted statistics** (Section 7.7) — weighted Pearson r, N_eff, p-values
9. **Orbital discriminator** (Section 8) — post-subtraction test on r̂ and k̂

This module is called by runner.py on each survey. It receives clean data (already quality-cut) and returns a complete results dict.

---

## TASK: BUILD boat/analysis.py

Create the file `boat/analysis.py` with the following content. Implement exactly as specified.

```python
"""
BOAT v2.0 — Statistical Analysis Engine

Implements the complete blind test protocol from the pre-registered manifest:
- Footprint correction (inverse local density weighting)
- 51-axis dipole regression (1 true + 50 decoys)
- Ratio statistic, rank test, p-value test
- Redshift-binned invariance test
- Dipole + quadrupole regression
- Orbital discriminator (post-subtraction residual test on r̂ and k̂)

Usage:
    from boat.analysis import run_full_analysis
    results = run_full_analysis(ra, dec, z_obs, survey_name)
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy import stats
import time

from boat.manifest import (
    V_OBS, C_KM_S, V_HAT, V_OBS_VEC, R_HAT, K_HAT,
    N_DECOY_AXES, DECOY_AXES_SEED, R_WEIGHT_DEG,
    P_VALUE_THRESHOLD, RANK_THRESHOLD,
    R_NULL_99, MAX_MIN_RATIO_ZBIN, MIN_BIN_SIZE,
    Z_BINS_STANDARD, Z_BINS_LRG,
    ORBITAL_DISC_P_THRESHOLD, ORBITAL_DISC_RATIO_THRESHOLD
)
from boat.geometry import radec_to_unit_vector


# ============================================================
# FOOTPRINT CORRECTION (Manifest Section 6)
# ============================================================

def compute_footprint_weights(ra, dec):
    """
    Compute inverse local density weights using cKDTree.
    
    For each galaxy:
    1. Count N_local = galaxies within R_weight (5°) angular radius
    2. Weight = 1 / N_local
    3. Normalise: weights / mean(weights)
    
    Uses chord distance on unit sphere for efficiency:
    chord = 2 * sin(θ/2) where θ = R_weight in radians
    
    Args:
        ra: array of RA in degrees
        dec: array of Dec in degrees
    
    Returns:
        weights: array of length N, normalised so mean = 1.0
    """
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    unit_vecs = np.column_stack([
        np.cos(dec_rad) * np.cos(ra_rad),
        np.cos(dec_rad) * np.sin(ra_rad),
        np.sin(dec_rad)
    ])
    
    chord_threshold = 2.0 * np.sin(np.radians(R_WEIGHT_DEG) / 2.0)
    
    tree = cKDTree(unit_vecs)
    counts = tree.query_ball_point(unit_vecs, r=chord_threshold, return_length=True)
    
    weights = 1.0 / counts.astype(np.float64)
    weights = weights / np.mean(weights)
    
    return weights


# ============================================================
# WEIGHTED STATISTICS (Manifest Section 7.7)
# ============================================================

def weighted_pearson_r(x, y, w):
    """
    Weighted Pearson correlation coefficient.
    
    r_w = Σ w_i (x_i - x̄_w)(y_i - ȳ_w) / 
          sqrt[Σ w_i (x_i - x̄_w)² × Σ w_i (y_i - ȳ_w)²]
    
    Args:
        x, y: data arrays
        w: weight array (mean ≈ 1)
    
    Returns:
        r: weighted correlation coefficient
    """
    w_sum = np.sum(w)
    x_mean = np.sum(w * x) / w_sum
    y_mean = np.sum(w * y) / w_sum
    
    dx = x - x_mean
    dy = y - y_mean
    
    numerator = np.sum(w * dx * dy)
    denom_x = np.sum(w * dx**2)
    denom_y = np.sum(w * dy**2)
    denominator = np.sqrt(denom_x * denom_y)
    
    if denominator == 0:
        return 0.0
    
    return numerator / denominator


def weighted_pearson_p(r, w):
    """
    P-value for weighted Pearson correlation.
    
    Uses effective sample size: N_eff = (Σ w_i)² / Σ w_i²
    Test statistic: t = r × sqrt(N_eff - 2) / sqrt(1 - r²)
    P-value from two-sided t-distribution with N_eff - 2 dof.
    
    Args:
        r: weighted correlation coefficient
        w: weight array
    
    Returns:
        p: two-sided p-value
    """
    n_eff = np.sum(w)**2 / np.sum(w**2)
    
    if n_eff <= 2 or abs(r) >= 1.0:
        return 1.0
    
    t_stat = r * np.sqrt(n_eff - 2) / np.sqrt(1 - r**2)
    p_value = 2.0 * stats.t.sf(abs(t_stat), df=n_eff - 2)
    
    return p_value


def weighted_linear_regression(x, y, w):
    """
    Weighted linear regression: y = A * x + intercept.
    
    Args:
        x, y: data arrays
        w: weight array
    
    Returns:
        (slope, intercept, r_weighted, p_value)
    """
    w_sum = np.sum(w)
    x_mean = np.sum(w * x) / w_sum
    y_mean = np.sum(w * y) / w_sum
    
    dx = x - x_mean
    dy = y - y_mean
    
    slope = np.sum(w * dx * dy) / np.sum(w * dx**2)
    intercept = y_mean - slope * x_mean
    
    r = weighted_pearson_r(x, y, w)
    p = weighted_pearson_p(r, w)
    
    return slope, intercept, r, p


# ============================================================
# DECOY AXIS GENERATION (Manifest Section 7.4)
# ============================================================

def generate_decoy_axes():
    """
    Generate 50 uniform random unit vectors on sphere.
    Seed: DECOY_AXES_SEED = 20260301 (locked in manifest).
    Identical to montecarlo.py implementation.
    
    Returns:
        axes: array of shape (50, 3)
    """
    rng = np.random.default_rng(DECOY_AXES_SEED)
    z = rng.uniform(-1, 1, size=N_DECOY_AXES)
    phi = rng.uniform(0, 2 * np.pi, size=N_DECOY_AXES)
    r_perp = np.sqrt(1 - z**2)
    
    return np.column_stack([
        r_perp * np.cos(phi),
        r_perp * np.sin(phi),
        z
    ])


# ============================================================
# UNIT VECTOR CONVERSION
# ============================================================

def ra_dec_to_unit_vectors(ra, dec):
    """Convert arrays of RA, Dec (degrees) to unit vectors (N, 3)."""
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    return np.column_stack([
        np.cos(dec_rad) * np.cos(ra_rad),
        np.cos(dec_rad) * np.sin(ra_rad),
        np.sin(dec_rad)
    ])


# ============================================================
# DIPOLE REGRESSION — PRIMARY TEST (Manifest Section 7.2)
# ============================================================

def run_dipole_test(dz, unit_vecs, weights, axes, axis_labels):
    """
    Run weighted dipole regression on all test axes.
    
    For each axis:
    1. Compute cos(θ) = unit_vec · axis
    2. Fit weighted linear model: Δz = A × cos(θ) + intercept
    3. Record slope (A), intercept, r, p-value
    
    Args:
        dz: residuals array (z_obs - mean(z_obs))
        unit_vecs: galaxy unit vectors (N, 3)
        weights: footprint weights (N,)
        axes: test axes (M, 3) — first row is true axis
        axis_labels: list of labels (M strings)
    
    Returns:
        list of dicts, one per axis, with keys:
        label, slope, intercept, r, abs_r, p_value, cos_theta
    """
    results = []
    cos_theta_all = axes @ unit_vecs.T  # shape (M, N)
    
    for i in range(len(axes)):
        cos_theta = cos_theta_all[i]
        slope, intercept, r, p = weighted_linear_regression(cos_theta, dz, weights)
        
        results.append({
            "label": axis_labels[i],
            "slope": slope,
            "intercept": intercept,
            "r": r,
            "abs_r": abs(r),
            "p_value": p,
            "cos_theta": cos_theta,  # kept for downstream use (orbital discriminator)
        })
    
    return results


# ============================================================
# DIPOLE + QUADRUPOLE REGRESSION (Manifest Section 7.3)
# ============================================================

def run_dipole_quadrupole(dz, cos_theta_true, weights):
    """
    Fit Δz = A × cos(θ) + B × (3cos²(θ) - 1)/2 + intercept
    on the true axis only.
    
    Args:
        dz: residuals
        cos_theta_true: cos(θ) for true axis
        weights: footprint weights
    
    Returns:
        dict with dipole_amplitude (A), quadrupole_amplitude (B),
        intercept, and p-values for each
    """
    # Design matrix: [cos(θ), (3cos²θ - 1)/2, 1]
    cos_t = cos_theta_true
    quad_term = (3.0 * cos_t**2 - 1.0) / 2.0
    ones = np.ones_like(cos_t)
    
    X = np.column_stack([cos_t, quad_term, ones])
    
    # Weighted least squares: (X^T W X)^-1 X^T W y
    W = np.diag(weights) if len(weights) < 50000 else None  # avoid huge matrix
    
    # For large N, use the normal equations directly
    XtW = X.T * weights[np.newaxis, :]  # (3, N) * (1, N) = (3, N)
    XtWX = XtW @ X  # (3, 3)
    XtWy = XtW @ dz  # (3,)
    
    try:
        coeffs = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return {
            "dipole_amplitude": 0.0,
            "quadrupole_amplitude": 0.0,
            "intercept": 0.0,
            "dipole_p": 1.0,
            "quadrupole_p": 1.0,
        }
    
    A, B, intercept = coeffs
    
    # P-values via residual analysis
    y_pred = X @ coeffs
    resid = dz - y_pred
    n_eff = np.sum(weights)**2 / np.sum(weights**2)
    
    # Weighted residual variance
    w_resid_var = np.sum(weights * resid**2) / (n_eff - 3)
    
    # Covariance matrix of coefficients
    try:
        cov = w_resid_var * np.linalg.inv(XtWX)
        se_A = np.sqrt(abs(cov[0, 0]))
        se_B = np.sqrt(abs(cov[1, 1]))
        
        t_A = A / se_A if se_A > 0 else 0.0
        t_B = B / se_B if se_B > 0 else 0.0
        
        p_A = 2.0 * stats.t.sf(abs(t_A), df=max(n_eff - 3, 1))
        p_B = 2.0 * stats.t.sf(abs(t_B), df=max(n_eff - 3, 1))
    except np.linalg.LinAlgError:
        p_A, p_B = 1.0, 1.0
    
    return {
        "dipole_amplitude": float(A),
        "quadrupole_amplitude": float(B),
        "intercept": float(intercept),
        "dipole_p": float(p_A),
        "quadrupole_p": float(p_B),
    }


# ============================================================
# REDSHIFT-BINNED INVARIANCE TEST (Manifest Section 7.6)
# ============================================================

def run_zbin_invariance(ra, dec, z_obs, weights, unit_vecs, survey_name):
    """
    Test that the dipole signal is consistent across redshift bins.
    
    Rules:
    - Use Z_BINS_STANDARD or Z_BINS_LRG depending on survey
    - Exclude bins with < MIN_BIN_SIZE galaxies
    - Need >= 2 qualifying bins, otherwise waive
    - True axis |r| must have same sign across qualifying bins
    - max(|r|) / min(|r|) <= MAX_MIN_RATIO_ZBIN
    
    Args:
        ra, dec, z_obs: arrays
        weights: footprint weights
        unit_vecs: galaxy unit vectors
        survey_name: string
    
    Returns:
        dict with pass/fail, per-bin results, waived flag
    """
    from boat.manifest import LRG_SURVEYS if hasattr(__import__('boat.manifest', fromlist=['LRG_SURVEYS']), 'LRG_SURVEYS') else None
    
    # Determine which bins to use
    lrg_surveys = {"DESI_LRG"}
    z_bins = Z_BINS_LRG if survey_name in lrg_surveys else Z_BINS_STANDARD
    
    true_axis = V_HAT
    cos_theta = unit_vecs @ true_axis
    
    bin_results = []
    qualifying_r_values = []
    
    for z_lo, z_hi in z_bins:
        mask = (z_obs >= z_lo) & (z_obs < z_hi)
        n_bin = np.sum(mask)
        
        if n_bin < MIN_BIN_SIZE:
            bin_results.append({
                "z_range": f"{z_lo}–{z_hi}",
                "n_galaxies": int(n_bin),
                "status": "excluded (< MIN_BIN_SIZE)",
                "r": None,
            })
            continue
        
        # Compute residuals and weighted r for this bin
        z_bin = z_obs[mask]
        dz_bin = z_bin - np.mean(z_bin)
        cos_bin = cos_theta[mask]
        w_bin = weights[mask]
        # Re-normalise bin weights
        w_bin = w_bin / np.mean(w_bin)
        
        r_bin = weighted_pearson_r(cos_bin, dz_bin, w_bin)
        p_bin = weighted_pearson_p(r_bin, w_bin)
        
        bin_results.append({
            "z_range": f"{z_lo}–{z_hi}",
            "n_galaxies": int(n_bin),
            "status": "qualifying",
            "r": float(r_bin),
            "abs_r": float(abs(r_bin)),
            "p_value": float(p_bin),
        })
        qualifying_r_values.append(r_bin)
    
    # Evaluate invariance
    n_qualifying = len(qualifying_r_values)
    
    if n_qualifying < 2:
        return {
            "pass": True,  # waived
            "waived": True,
            "reason": f"Only {n_qualifying} qualifying bins (need >= 2)",
            "n_qualifying_bins": n_qualifying,
            "bins": bin_results,
        }
    
    # Check same sign
    signs = [np.sign(r) for r in qualifying_r_values]
    same_sign = all(s == signs[0] for s in signs)
    
    # Check max/min ratio
    abs_r_vals = [abs(r) for r in qualifying_r_values]
    if min(abs_r_vals) > 0:
        ratio = max(abs_r_vals) / min(abs_r_vals)
    else:
        ratio = float('inf')
    
    invariance_pass = same_sign and (ratio <= MAX_MIN_RATIO_ZBIN)
    
    return {
        "pass": invariance_pass,
        "waived": False,
        "same_sign": same_sign,
        "max_min_ratio": round(float(ratio), 4),
        "threshold": MAX_MIN_RATIO_ZBIN,
        "n_qualifying_bins": n_qualifying,
        "bins": bin_results,
    }


# ============================================================
# ORBITAL DISCRIMINATOR (Manifest Section 8)
# ============================================================

def run_orbital_discriminator(dz, unit_vecs, weights, primary_result):
    """
    After subtracting the primary v̂ dipole, test residuals against
    r̂ (CoG X direction) and k̂ (orbital normal).
    
    Steps:
    1. Subtract fitted v̂ dipole: Δz_resid = Δz - A_v × cos(θ_v) - intercept
    2. Compute weighted r and p against r̂ and k̂
    3. Compare against 50 decoy axes (same post-subtraction test)
    4. Report DETECTION or NON-DETECTION
    
    Args:
        dz: original residuals (z_obs - mean(z_obs))
        unit_vecs: galaxy unit vectors (N, 3)
        weights: footprint weights
        primary_result: dict from run_dipole_test for the true axis
            (must contain 'slope', 'intercept', 'cos_theta')
    
    Returns:
        dict with r̂ and k̂ test results, detection verdict
    """
    # Step 1: Subtract primary dipole
    cos_theta_v = primary_result["cos_theta"]  # cos(θ) w.r.t. v̂
    A_v = primary_result["slope"]
    intercept_v = primary_result["intercept"]
    
    dz_residual = dz - A_v * cos_theta_v - intercept_v
    
    # Step 2: Test against r̂ and k̂
    test_axes = {"r_hat": R_HAT, "k_hat": K_HAT}
    discriminator_results = {}
    
    for axis_name, axis_vec in test_axes.items():
        cos_theta = unit_vecs @ axis_vec
        slope, intercept, r, p = weighted_linear_regression(cos_theta, dz_residual, weights)
        
        discriminator_results[axis_name] = {
            "r": float(r),
            "abs_r": float(abs(r)),
            "p_value": float(p),
            "slope": float(slope),
            "intercept": float(intercept),
        }
    
    # Step 3: Compare against 50 decoy axes (post-subtraction)
    decoy_axes = generate_decoy_axes()
    cos_theta_decoys = decoy_axes @ unit_vecs.T  # (50, N)
    
    decoy_abs_r = []
    for j in range(N_DECOY_AXES):
        r_decoy = weighted_pearson_r(cos_theta_decoys[j], dz_residual, weights)
        decoy_abs_r.append(abs(r_decoy))
    
    mean_decoy_abs_r = np.mean(decoy_abs_r)
    
    # Compute ratios for r̂ and k̂
    for axis_name in ["r_hat", "k_hat"]:
        abs_r = discriminator_results[axis_name]["abs_r"]
        if mean_decoy_abs_r > 0:
            ratio = abs_r / mean_decoy_abs_r
        else:
            ratio = 0.0
        discriminator_results[axis_name]["ratio"] = float(ratio)
    
    # Step 4: Evaluate DETECTION
    detection = False
    detected_axes = []
    
    for axis_name in ["r_hat", "k_hat"]:
        res = discriminator_results[axis_name]
        p_pass = res["p_value"] < ORBITAL_DISC_P_THRESHOLD
        ratio_pass = res["ratio"] >= ORBITAL_DISC_RATIO_THRESHOLD
        res["p_pass"] = p_pass
        res["ratio_pass"] = ratio_pass
        res["detected"] = p_pass and ratio_pass
        if res["detected"]:
            detection = True
            detected_axes.append(axis_name)
    
    discriminator_results["mean_decoy_abs_r"] = float(mean_decoy_abs_r)
    discriminator_results["detection"] = detection
    discriminator_results["detected_axes"] = detected_axes
    discriminator_results["verdict"] = (
        f"DETECTED on {', '.join(detected_axes)}" if detection
        else "NOT DETECTED"
    )
    
    return discriminator_results


# ============================================================
# z_model COMPUTATION (for reference only — NOT in pipeline)
# ============================================================

def compute_z_model(ra, dec):
    """
    Compute predicted redshift from orbital motion for each galaxy.
    
    1. Convert (RA, Dec) → unit vector n̂
    2. v_LOS = -V_OBS_VEC · n̂
    3. β = v_LOS / c
    4. z_model = sqrt((1+β)/(1-β)) - 1
    
    THIS IS FOR REPORTING ONLY. z_model never enters the statistical
    analysis pipeline. It is compared against the observed dipole as
    a consistency check.
    
    Args:
        ra, dec: arrays in degrees
    
    Returns:
        z_model: array of predicted redshifts
    """
    unit_vecs = ra_dec_to_unit_vectors(ra, dec)
    v_los = -unit_vecs @ V_OBS_VEC  # (N,) dot product with velocity vector
    beta = v_los / C_KM_S
    z_model = np.sqrt((1 + beta) / (1 - beta)) - 1
    return z_model


# ============================================================
# MAIN ANALYSIS FUNCTION
# ============================================================

def run_full_analysis(ra, dec, z_obs, survey_name):
    """
    Run the complete BOAT v2.0 analysis protocol on one survey.
    
    Args:
        ra: array of RA in degrees (after quality cuts)
        dec: array of Dec in degrees
        z_obs: array of observed redshifts
        survey_name: string identifying the survey
    
    Returns:
        dict with complete results including:
        - primary dipole test (all 51 axes)
        - ratio, rank, p-value verdicts
        - z-bin invariance test
        - dipole + quadrupole regression
        - orbital discriminator
        - z_model comparison (reference only)
        - per-dataset PASS/FAIL verdict
    """
    t_start = time.time()
    n_galaxies = len(ra)
    
    print(f"\n{'='*60}")
    print(f"BOAT v2.0 Analysis: {survey_name}")
    print(f"N = {n_galaxies:,} galaxies")
    print(f"{'='*60}")
    
    # ---- Footprint weights ----
    print(f"\n  Computing footprint weights (R = {R_WEIGHT_DEG}°)...")
    t0 = time.time()
    weights = compute_footprint_weights(ra, dec)
    print(f"  Weights done in {time.time()-t0:.1f}s")
    
    # ---- Unit vectors ----
    unit_vecs = ra_dec_to_unit_vectors(ra, dec)
    
    # ---- Residuals ----
    dz = z_obs - np.mean(z_obs)
    
    # ---- Generate axes ----
    decoy_axes = generate_decoy_axes()
    true_axis = V_HAT.reshape(1, 3)
    all_axes = np.vstack([true_axis, decoy_axes])
    axis_labels = ["TRUE (CMB dipole)"] + [f"DECOY_{i+1:02d}" for i in range(N_DECOY_AXES)]
    
    # ---- Primary dipole test on all 51 axes ----
    print(f"  Running dipole regression on {len(all_axes)} axes...")
    dipole_results = run_dipole_test(dz, unit_vecs, weights, all_axes, axis_labels)
    
    # Extract true axis result
    true_result = dipole_results[0]
    decoy_results = dipole_results[1:]
    
    # ---- Ratio statistic ----
    decoy_abs_r = [d["abs_r"] for d in decoy_results]
    mean_decoy_abs_r = np.mean(decoy_abs_r)
    
    if mean_decoy_abs_r > 0:
        ratio = true_result["abs_r"] / mean_decoy_abs_r
    else:
        ratio = 0.0
    
    # ---- Rank ----
    all_abs_r = [d["abs_r"] for d in dipole_results]
    sorted_indices = np.argsort(all_abs_r)[::-1]
    rank = int(np.where(sorted_indices == 0)[0][0]) + 1  # 1-indexed
    
    # ---- PASS criteria ----
    ratio_pass = (R_NULL_99 is not None) and (ratio >= R_NULL_99)
    rank_pass = rank <= RANK_THRESHOLD
    p_pass = true_result["p_value"] < P_VALUE_THRESHOLD
    
    print(f"\n  True axis results:")
    print(f"    |r| = {true_result['abs_r']:.6f}")
    print(f"    p-value = {true_result['p_value']:.2e}")
    print(f"    Ratio = {ratio:.4f} (threshold: {R_NULL_99})")
    print(f"    Rank = #{rank} of {len(all_axes)} (threshold: top {RANK_THRESHOLD})")
    print(f"    Ratio PASS: {ratio_pass}")
    print(f"    Rank PASS: {rank_pass}")
    print(f"    p-value PASS: {p_pass}")
    
    # ---- z-bin invariance ----
    print(f"\n  Running z-bin invariance test...")
    zbin_result = run_zbin_invariance(ra, dec, z_obs, weights, unit_vecs, survey_name)
    zbin_pass = zbin_result["pass"]
    print(f"    z-bin invariance PASS: {zbin_pass}")
    if not zbin_result["waived"]:
        print(f"    Same sign: {zbin_result.get('same_sign')}")
        print(f"    Max/min ratio: {zbin_result.get('max_min_ratio')}")
    
    # ---- Overall per-dataset verdict ----
    if R_NULL_99 is None:
        overall_pass = False
        pass_note = "R_NULL_99 not yet set — cannot evaluate"
    else:
        overall_pass = ratio_pass and rank_pass and p_pass and zbin_pass
    
    print(f"\n  *** {survey_name} VERDICT: {'PASS' if overall_pass else 'FAIL'} ***")
    
    # ---- Dipole + quadrupole regression ----
    print(f"\n  Running dipole + quadrupole regression...")
    dq_result = run_dipole_quadrupole(dz, true_result["cos_theta"], weights)
    print(f"    Dipole amplitude (A): {dq_result['dipole_amplitude']:.6e}")
    print(f"    Quadrupole amplitude (B): {dq_result['quadrupole_amplitude']:.6e}")
    
    # ---- Orbital discriminator ----
    print(f"\n  Running orbital discriminator (post-subtraction test)...")
    orbital_result = run_orbital_discriminator(dz, unit_vecs, weights, true_result)
    print(f"    r̂ (CoG X): r = {orbital_result['r_hat']['r']:.6f}, "
          f"p = {orbital_result['r_hat']['p_value']:.2e}, "
          f"ratio = {orbital_result['r_hat']['ratio']:.4f}")
    print(f"    k̂ (normal): r = {orbital_result['k_hat']['r']:.6f}, "
          f"p = {orbital_result['k_hat']['p_value']:.2e}, "
          f"ratio = {orbital_result['k_hat']['ratio']:.4f}")
    print(f"    Orbital discriminator: {orbital_result['verdict']}")
    
    # ---- z_model reference ----
    z_model = compute_z_model(ra, dec)
    
    elapsed = time.time() - t_start
    print(f"\n  Analysis complete in {elapsed:.1f}s")
    
    # ---- Compile results ----
    # Strip cos_theta from dipole results (large arrays, not needed in output)
    dipole_output = []
    for d in dipole_results:
        d_out = {k: v for k, v in d.items() if k != "cos_theta"}
        dipole_output.append(d_out)
    
    results = {
        "survey": survey_name,
        "n_galaxies": n_galaxies,
        "elapsed_seconds": round(elapsed, 1),
        "mean_z_obs": float(np.mean(z_obs)),
        "mean_z_model_ref": float(np.mean(z_model)),
        "weight_stats": {
            "min": float(weights.min()),
            "max": float(weights.max()),
            "mean": float(weights.mean()),
            "std": float(weights.std()),
        },
        "primary_test": {
            "true_axis": {
                "r": float(true_result["r"]),
                "abs_r": float(true_result["abs_r"]),
                "p_value": float(true_result["p_value"]),
                "slope": float(true_result["slope"]),
                "intercept": float(true_result["intercept"]),
            },
            "ratio": float(ratio),
            "rank": rank,
            "mean_decoy_abs_r": float(mean_decoy_abs_r),
            "R_NULL_99": R_NULL_99,
            "ratio_pass": ratio_pass,
            "rank_pass": rank_pass,
            "p_pass": p_pass,
        },
        "zbin_invariance": zbin_result,
        "dipole_quadrupole": dq_result,
        "orbital_discriminator": orbital_result,
        "verdict": "PASS" if overall_pass else "FAIL",
        "verdict_note": "" if R_NULL_99 is not None else "R_NULL_99 not set",
        "all_axes": dipole_output,
    }
    
    return results
```

---

## IMPLEMENTATION NOTES

### 1. The z-bin invariance function has a conditional import

The line checking for `LRG_SURVEYS` in manifest is awkward. Replace it with a simple hardcoded set:

```python
lrg_surveys = {"DESI_LRG"}
```

This matches what's in montecarlo.py. The set is small and unlikely to change.

### 2. cos_theta is kept in memory during analysis

The `run_dipole_test` function stores `cos_theta` in each result dict for use by the orbital discriminator. This is stripped out before the final output (the compile step removes it). This avoids recomputing the dot products.

### 3. R_NULL_99 = None handling

When R_NULL_99 is still `None` (before sealing), the ratio test automatically fails and the verdict shows a note. This lets us test the full pipeline during development without a sealed manifest.

### 4. z_model is computed but never used in statistics

`compute_z_model` produces predicted redshifts for reference/reporting only. It appears in the output as `mean_z_model_ref`. It is NEVER subtracted from data, NEVER used in any correlation, NEVER enters the pipeline. This is a v1.1 lesson (Critical Error #1).

---

## VERIFICATION TESTS

After building, run these:

### Test 1: Quick analysis on 6dFGS
```python
from boat.montecarlo import load_6dfgs, apply_quality_cuts
from boat.analysis import run_full_analysis

df = load_6dfgs()
df = apply_quality_cuts(df, "6dFGS")

results = run_full_analysis(
    df['ra'].values, df['dec'].values, df['z'].values, "6dFGS"
)

print(f"\n6dFGS verdict: {results['verdict']}")
print(f"True axis |r|: {results['primary_test']['true_axis']['abs_r']:.6f}")
print(f"Ratio: {results['primary_test']['ratio']:.4f}")
print(f"Rank: #{results['primary_test']['rank']}")
print(f"Orbital discriminator: {results['orbital_discriminator']['verdict']}")
```

This should complete in under 30 seconds. It will show FAIL for the verdict because R_NULL_99 is None, but all the numbers should be populated.

### Test 2: Verify v1.1 consistency

6dFGS was the one survey that PASSED v1.1. Its |r| was 0.077 and it ranked #1 of 11 axes. In v2 (with 51 axes and footprint weighting), the |r| and rank may differ, but we should see the true axis performing strongly.

**Expected behaviour:**
- |r| should be in a similar range (order of magnitude 0.01–0.1)
- True axis should rank near the top
- p-value should be very small

---

## REPORT BACK

1. analysis.py created: yes/no
2. Test 1 (6dFGS quick analysis):
   - Completed: yes/no
   - Time taken
   - True axis |r|, ratio, rank, p-value
   - z-bin invariance result
   - Orbital discriminator verdict
   - Overall verdict (expected: FAIL due to R_NULL_99 = None)
3. Any errors or warnings

**Do not run on SDSS yet.** Wait for instruction.

---

*Build prompt prepared by Claude Opus (Anthropic), 22 February 2026*
