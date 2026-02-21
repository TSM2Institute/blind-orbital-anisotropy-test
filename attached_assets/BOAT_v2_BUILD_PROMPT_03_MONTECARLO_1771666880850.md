# BOAT v2.0 — Build Prompt #3: montecarlo.py

**Context:** manifest.py is built and verified (all 8 checks pass). geometry.py imports cleanly. Data files renamed. We now build the Monte Carlo null calibration module — the most computationally intensive part of the pipeline.

---

## WHAT THIS MODULE DOES

montecarlo.py determines R_NULL_99: the ratio threshold that ensures a global false-positive rate < 0.01 under the null hypothesis (no real dipole signal).

**Method per survey:**
1. Load real survey data (positions + redshifts)
2. Apply quality cuts (z range, galactic plane mask, duplicates)
3. Compute footprint weights ONCE (positions are fixed)
4. Precompute cos(θ) arrays for all 51 axes ONCE
5. For each of 10,000 realisations:
   a. Randomly shuffle z_obs among positions (destroys any real dipole)
   b. Compute Δz = z_shuffled − mean(z_shuffled)
   c. Compute weighted Pearson r for all 51 axes
   d. Record ratio = |r_true| / mean(|r_decoys|)
6. Combine ratio values across all surveys
7. R_NULL_99 = 99th percentile of combined distribution

---

## ARCHITECTURE NOTES

### Performance-critical precomputations

The key insight: since we keep (RA, Dec) positions fixed and only shuffle z values, three expensive computations happen ONCE per survey, not per realisation:

1. **Footprint weights** — depend only on positions
2. **cos(θ) arrays** — dot product of each galaxy unit vector with each axis, depend only on positions
3. **Decoy axes** — generated from fixed seed

The inner loop (per realisation) then reduces to: shuffle an array, subtract its mean, compute 51 weighted correlations using precomputed arrays. This is fast.

### Footprint weighting

For each galaxy i, count how many other galaxies fall within R_weight = 5° angular radius. Weight = 1/count, then normalise so mean(weights) = 1.

**For large surveys (SDSS: ~1.5M galaxies), brute-force pairwise angular distance is O(N²) and will not complete.** Use scipy.spatial.cKDTree on unit vectors with a Euclidean distance threshold equivalent to 5° angular separation:

- Angular separation θ = 5°
- Equivalent chord distance: d = 2 × sin(θ/2) = 2 × sin(2.5°) ≈ 0.08716

This converts the O(N²) problem to approximately O(N log N).

### DESI surveys

We don't have DESI data yet. The module must support:
- Running on 6dFGS + SDSS now (2 surveys)
- Adding DESI BGS + LRG later (after download)
- Combining all results before sealing

Design the module so it processes whatever surveys are available in `attached_assets/` and stores per-survey results that can be merged later.

---

## TASK: BUILD boat/montecarlo.py

Create the file `boat/montecarlo.py` with the following structure. I am providing the architecture and key algorithms. Implement them exactly.

```python
"""
BOAT v2.0 — Monte Carlo Null Calibration

Determines R_NULL_99 by running the complete 51-axis protocol on
shuffled (null) versions of each survey. Footprint weights and
cos(θ) arrays are precomputed once per survey for performance.

Usage:
    from boat.montecarlo import run_null_calibration
    results = run_null_calibration(surveys=["6dFGS", "SDSS"])
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy import stats
import json
import time
import os

from boat.manifest import (
    V_OBS, C_KM_S, V_HAT, R_HAT, K_HAT, V_OBS_VEC,
    N_DECOY_AXES, DECOY_AXES_SEED, R_WEIGHT_DEG,
    Z_MIN, Z_MAX, Z_MIN_LRG, Z_MAX_LRG,
    GALACTIC_PLANE_MASK_DEG, MIN_SAMPLE_SIZE,
    N_MC_REALISATIONS, DATA_DIR, DATA_FILES,
    MIN_SKY_FRACTION, MIN_MEDIAN_Z
)


# ============================================================
# SURVEY LOADING (minimal — just enough for Monte Carlo)
# ============================================================

def load_6dfgs():
    """
    Load 6dFGS DR3 from gzipped text file.
    
    File format: space-separated, columns include RA (RAJ2000), 
    Dec (DEJ2000), and redshift (cz converted to z, or zfinal).
    
    You need to inspect the actual file header to determine exact
    column names. The file is at: attached_assets/6dFGS_DR3.txt.gz
    
    Read with pandas, handling the text format appropriately.
    Return DataFrame with columns: 'ra', 'dec', 'z'
    """
    import pandas as pd
    filepath = os.path.join(DATA_DIR, DATA_FILES["6dFGS"])
    
    # 6dFGS DR3 format: pipe-delimited (|) with header rows
    # Key columns: RAJ2000, DEJ2000, zfinal
    # IMPORTANT: Inspect the actual file to confirm column names
    # and delimiter before finalising this loader.
    #
    # Read the first few lines of the file to determine format:
    import gzip
    with gzip.open(filepath, 'rt') as f:
        for i in range(5):
            print(f"Line {i}: {f.readline().strip()}")
    
    # TODO: Complete this loader based on actual file inspection.
    # Must return DataFrame with columns: 'ra', 'dec', 'z'
    raise NotImplementedError("Inspect file format and complete loader")


def load_sdss():
    """
    Load SDSS DR18 from CSV file.
    
    File is a CasJobs export. Expected columns include:
    ra, dec, z (or similar names).
    
    File at: attached_assets/SDSS_DR18.csv
    
    Read with pandas. Return DataFrame with columns: 'ra', 'dec', 'z'
    """
    import pandas as pd
    filepath = os.path.join(DATA_DIR, DATA_FILES["SDSS"])
    
    # Read first few lines to confirm column names
    df_head = pd.read_csv(filepath, nrows=3)
    print(f"SDSS columns: {list(df_head.columns)}")
    print(df_head.head())
    
    # TODO: Complete this loader based on actual column names.
    # Must return DataFrame with columns: 'ra', 'dec', 'z'
    raise NotImplementedError("Inspect column names and complete loader")


SURVEY_LOADERS = {
    "6dFGS": load_6dfgs,
    "SDSS": load_sdss,
    # "DESI_BGS": load_desi_bgs,  # Added after download
    # "DESI_LRG": load_desi_lrg,  # Added after download
}

# Identify which surveys use LRG z-range
LRG_SURVEYS = {"DESI_LRG"}


# ============================================================
# QUALITY CUTS
# ============================================================

def apply_quality_cuts(df, survey_name):
    """
    Apply standard quality cuts to a survey DataFrame.
    
    Args:
        df: DataFrame with columns 'ra', 'dec', 'z'
        survey_name: string identifying the survey
    
    Returns:
        DataFrame after cuts, with report printed
    """
    n_start = len(df)
    
    # 1. Redshift range
    if survey_name in LRG_SURVEYS:
        z_min, z_max = Z_MIN_LRG, Z_MAX_LRG
    else:
        z_min, z_max = Z_MIN, Z_MAX
    
    df = df[(df['z'] >= z_min) & (df['z'] <= z_max)].copy()
    n_after_z = len(df)
    
    # 2. Galactic plane mask (|b| >= 10°)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coords = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg, frame='icrs')
    galactic = coords.galactic
    mask_gal = np.abs(galactic.b.deg) >= GALACTIC_PLANE_MASK_DEG
    df = df[mask_gal].copy()
    n_after_gal = len(df)
    
    # 3. Duplicate removal (nearest neighbour < 1 arcsec)
    from astropy.coordinates import match_coordinates_sky
    coords_clean = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg, frame='icrs')
    idx, sep, _ = match_coordinates_sky(coords_clean, coords_clean, nthneighbor=2)
    mask_dup = sep.arcsec >= 1.0
    df = df[mask_dup].copy()
    n_after_dup = len(df)
    
    print(f"  Quality cuts for {survey_name}:")
    print(f"    Start: {n_start:,}")
    print(f"    After z-cut ({z_min}–{z_max}): {n_after_z:,}")
    print(f"    After galactic mask (|b| >= {GALACTIC_PLANE_MASK_DEG}°): {n_after_gal:,}")
    print(f"    After duplicate removal: {n_after_dup:,}")
    
    return df


# ============================================================
# PRE-QUALIFICATION CHECK
# ============================================================

def check_prequalification(df, survey_name):
    """
    Check if survey meets pre-qualification criteria.
    
    Returns:
        (passes: bool, sky_fraction: float, median_z: float)
    """
    # Sky fraction: count occupied 1° × 1° bins
    ra_bins = np.floor(df['ra'].values).astype(int)  # 0–359
    dec_bins = np.floor(df['dec'].values + 90).astype(int)  # 0–179
    bin_ids = ra_bins * 180 + dec_bins
    n_occupied = len(np.unique(bin_ids))
    total_bins = 360 * 180  # 64,800
    sky_fraction = n_occupied / total_bins
    
    median_z = np.median(df['z'].values)
    
    passes = (sky_fraction >= MIN_SKY_FRACTION) and (median_z >= MIN_MEDIAN_Z) and (len(df) >= MIN_SAMPLE_SIZE)
    
    print(f"  Pre-qualification for {survey_name}:")
    print(f"    Sky fraction: {sky_fraction:.1%} (threshold: {MIN_SKY_FRACTION:.0%}) — {'PASS' if sky_fraction >= MIN_SKY_FRACTION else 'FAIL'}")
    print(f"    Median z: {median_z:.4f} (threshold: {MIN_MEDIAN_Z}) — {'PASS' if median_z >= MIN_MEDIAN_Z else 'FAIL'}")
    print(f"    Sample size: {len(df):,} (threshold: {MIN_SAMPLE_SIZE:,}) — {'PASS' if len(df) >= MIN_SAMPLE_SIZE else 'FAIL'}")
    print(f"    Overall: {'QUALIFIED' if passes else 'DISQUALIFIED'}")
    
    return passes, sky_fraction, median_z


# ============================================================
# FOOTPRINT WEIGHTING
# ============================================================

def compute_footprint_weights(ra, dec):
    """
    Compute inverse local density weights.
    
    For each galaxy, count neighbours within R_weight (5°) using
    cKDTree on unit vectors with chord distance threshold.
    
    Args:
        ra: array of RA in degrees
        dec: array of Dec in degrees
    
    Returns:
        weights: array, normalised so mean = 1.0
    """
    # Convert to unit vectors
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    unit_vecs = np.column_stack([x, y, z])
    
    # Chord distance threshold for R_weight degrees
    chord_threshold = 2.0 * np.sin(np.radians(R_WEIGHT_DEG) / 2.0)
    
    # Build KD-tree and count neighbours
    tree = cKDTree(unit_vecs)
    counts = tree.query_ball_point(unit_vecs, r=chord_threshold, return_length=True)
    
    # counts includes the point itself, so minimum is 1
    # Weight = 1 / count, normalised to mean = 1
    weights = 1.0 / counts.astype(np.float64)
    weights = weights / np.mean(weights)
    
    return weights


# ============================================================
# DECOY AXIS GENERATION
# ============================================================

def generate_decoy_axes():
    """
    Generate N_DECOY_AXES uniform random unit vectors on sphere.
    
    Uses the locked seed DECOY_AXES_SEED = 20260301.
    Same axes used for all surveys and all realisations.
    
    Returns:
        decoy_axes: array of shape (N_DECOY_AXES, 3)
    """
    rng = np.random.default_rng(DECOY_AXES_SEED)
    
    # Uniform on sphere: z uniform in [-1, 1], phi uniform in [0, 2π)
    z = rng.uniform(-1, 1, size=N_DECOY_AXES)
    phi = rng.uniform(0, 2 * np.pi, size=N_DECOY_AXES)
    r_perp = np.sqrt(1 - z**2)
    
    axes = np.column_stack([
        r_perp * np.cos(phi),
        r_perp * np.sin(phi),
        z
    ])
    
    return axes


# ============================================================
# WEIGHTED PEARSON CORRELATION
# ============================================================

def weighted_pearson_r(x, y, w):
    """
    Compute weighted Pearson correlation coefficient.
    
    r_w = Σ w_i (x_i - x̄_w)(y_i - ȳ_w) / 
          sqrt(Σ w_i (x_i - x̄_w)² × Σ w_i (y_i - ȳ_w)²)
    
    Args:
        x, y: arrays of same length
        w: weight array (should have mean ≈ 1)
    
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


# ============================================================
# CORE MONTE CARLO LOOP
# ============================================================

def run_mc_for_survey(ra, dec, z_obs, weights, cos_theta_true, cos_theta_decoys, n_realisations, rng_seed):
    """
    Run Monte Carlo null simulation for one survey.
    
    All position-dependent quantities (weights, cos_theta arrays) are
    precomputed and passed in. The only per-realisation operation is
    shuffling z_obs and computing weighted correlations.
    
    Args:
        ra: array of RA (not used in loop — included for API clarity)
        dec: array of Dec (not used in loop)
        z_obs: array of observed redshifts
        weights: footprint weights (precomputed, mean = 1)
        cos_theta_true: array of cos(θ) for true axis, shape (N,)
        cos_theta_decoys: array of cos(θ) for decoy axes, shape (N_DECOY, N)
        n_realisations: number of MC realisations
        rng_seed: seed for this survey's shuffling RNG
    
    Returns:
        ratios: array of shape (n_realisations,) — |r_true|/mean(|r_decoy|) per realisation
    """
    rng = np.random.default_rng(rng_seed)
    n_galaxies = len(z_obs)
    n_decoys = cos_theta_decoys.shape[0]
    ratios = np.empty(n_realisations)
    
    for i in range(n_realisations):
        # Shuffle redshifts (destroys any real dipole)
        z_shuffled = rng.permutation(z_obs)
        
        # Compute residuals
        dz = z_shuffled - np.mean(z_shuffled)
        
        # Weighted correlation with true axis
        r_true = weighted_pearson_r(cos_theta_true, dz, weights)
        
        # Weighted correlation with each decoy axis
        r_decoys = np.empty(n_decoys)
        for j in range(n_decoys):
            r_decoys[j] = weighted_pearson_r(cos_theta_decoys[j], dz, weights)
        
        # Ratio statistic
        mean_decoy_abs_r = np.mean(np.abs(r_decoys))
        if mean_decoy_abs_r == 0:
            ratios[i] = 0.0
        else:
            ratios[i] = np.abs(r_true) / mean_decoy_abs_r
        
        # Progress reporting every 1000 realisations
        if (i + 1) % 1000 == 0:
            print(f"    Realisation {i+1:,}/{n_realisations:,} — "
                  f"latest ratio: {ratios[i]:.4f}")
    
    return ratios


# ============================================================
# PRECOMPUTATION
# ============================================================

def precompute_cos_theta(unit_vectors, axes):
    """
    Compute cos(θ) between each galaxy and each axis.
    
    Args:
        unit_vectors: galaxy unit vectors, shape (N, 3)
        axes: test axes, shape (M, 3)
    
    Returns:
        cos_theta: array of shape (M, N)
    """
    # Matrix multiplication: (M, 3) @ (3, N) = (M, N)
    return axes @ unit_vectors.T


def ra_dec_to_unit_vectors(ra, dec):
    """Convert arrays of RA, Dec (degrees) to unit vectors, shape (N, 3)."""
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    return np.column_stack([
        np.cos(dec_rad) * np.cos(ra_rad),
        np.cos(dec_rad) * np.sin(ra_rad),
        np.sin(dec_rad)
    ])


# ============================================================
# MAIN CALIBRATION FUNCTION
# ============================================================

def run_null_calibration(surveys=None, n_realisations=None, mc_base_seed=42):
    """
    Run the full Monte Carlo null calibration.
    
    Args:
        surveys: list of survey names to process (default: all available)
        n_realisations: override N_MC_REALISATIONS (useful for testing)
        mc_base_seed: base seed for MC shuffling RNG (each survey gets
                      mc_base_seed + survey_index)
    
    Returns:
        dict with complete results including R_NULL_99 and per-survey data
    """
    if n_realisations is None:
        n_realisations = N_MC_REALISATIONS
    
    if surveys is None:
        surveys = list(SURVEY_LOADERS.keys())
    
    print("=" * 60)
    print("BOAT v2.0 — Monte Carlo Null Calibration")
    print(f"Surveys: {surveys}")
    print(f"Realisations per survey: {n_realisations:,}")
    print("=" * 60)
    
    # Generate decoy axes (fixed for all surveys)
    decoy_axes = generate_decoy_axes()
    true_axis = V_HAT.reshape(1, 3)  # CMB dipole direction
    all_axes = np.vstack([true_axis, decoy_axes])  # shape (51, 3)
    
    print(f"\nDecoy axes generated: {len(decoy_axes)} (seed: {DECOY_AXES_SEED})")
    print(f"Total axes per dataset: {len(all_axes)}")
    
    # Process each survey
    all_ratios = []
    per_survey_results = {}
    
    for s_idx, survey_name in enumerate(surveys):
        print(f"\n{'='*60}")
        print(f"Processing: {survey_name}")
        print(f"{'='*60}")
        
        t_start = time.time()
        
        # Load data
        loader = SURVEY_LOADERS[survey_name]
        df = loader()
        print(f"  Loaded: {len(df):,} galaxies")
        
        # Quality cuts
        df = apply_quality_cuts(df, survey_name)
        
        # Pre-qualification check
        passes, sky_frac, med_z = check_prequalification(df, survey_name)
        if not passes:
            print(f"  SKIPPING {survey_name} — does not pre-qualify")
            per_survey_results[survey_name] = {
                "status": "disqualified",
                "sky_fraction": sky_frac,
                "median_z": med_z,
                "n_galaxies": len(df)
            }
            continue
        
        ra = df['ra'].values
        dec = df['dec'].values
        z_obs = df['z'].values
        
        # Precompute footprint weights (ONCE)
        print(f"\n  Computing footprint weights (R = {R_WEIGHT_DEG}°)...")
        t_weights = time.time()
        weights = compute_footprint_weights(ra, dec)
        print(f"  Weights computed in {time.time() - t_weights:.1f}s")
        print(f"  Weight stats: min={weights.min():.4f}, max={weights.max():.4f}, "
              f"mean={weights.mean():.4f}, std={weights.std():.4f}")
        
        # Precompute unit vectors and cos(θ) arrays (ONCE)
        unit_vecs = ra_dec_to_unit_vectors(ra, dec)
        cos_theta_all = precompute_cos_theta(unit_vecs, all_axes)  # shape (51, N)
        cos_theta_true = cos_theta_all[0]       # shape (N,)
        cos_theta_decoys = cos_theta_all[1:]    # shape (50, N)
        
        # Run MC
        print(f"\n  Running {n_realisations:,} null realisations...")
        survey_seed = mc_base_seed + s_idx
        ratios = run_mc_for_survey(
            ra, dec, z_obs, weights,
            cos_theta_true, cos_theta_decoys,
            n_realisations, survey_seed
        )
        
        elapsed = time.time() - t_start
        
        # Per-survey statistics
        p50 = np.percentile(ratios, 50)
        p95 = np.percentile(ratios, 95)
        p99 = np.percentile(ratios, 99)
        p995 = np.percentile(ratios, 99.5)
        
        print(f"\n  {survey_name} complete in {elapsed:.1f}s")
        print(f"  Null ratio distribution:")
        print(f"    Median (50th): {p50:.4f}")
        print(f"    95th percentile: {p95:.4f}")
        print(f"    99th percentile: {p99:.4f}")
        print(f"    99.5th percentile: {p995:.4f}")
        print(f"    Min: {ratios.min():.4f}, Max: {ratios.max():.4f}")
        
        all_ratios.extend(ratios.tolist())
        
        per_survey_results[survey_name] = {
            "status": "completed",
            "n_galaxies": len(df),
            "sky_fraction": sky_frac,
            "median_z": med_z,
            "n_realisations": n_realisations,
            "seed": survey_seed,
            "elapsed_seconds": round(elapsed, 1),
            "percentiles": {
                "p50": round(p50, 6),
                "p95": round(p95, 6),
                "p99": round(p99, 6),
                "p995": round(p995, 6),
            },
            "min_ratio": round(float(ratios.min()), 6),
            "max_ratio": round(float(ratios.max()), 6),
            "mean_ratio": round(float(ratios.mean()), 6),
            "weight_stats": {
                "min": round(float(weights.min()), 6),
                "max": round(float(weights.max()), 6),
                "mean": round(float(weights.mean()), 6),
                "std": round(float(weights.std()), 6),
            }
        }
    
    # ============================================================
    # COMBINED RESULTS
    # ============================================================
    
    all_ratios = np.array(all_ratios)
    n_total = len(all_ratios)
    
    if n_total == 0:
        print("\nERROR: No surveys completed. Cannot compute R_NULL_99.")
        return {"error": "No surveys completed"}
    
    R_NULL_95 = float(np.percentile(all_ratios, 95))
    R_NULL_99 = float(np.percentile(all_ratios, 99))
    R_NULL_995 = float(np.percentile(all_ratios, 99.5))
    
    print(f"\n{'='*60}")
    print(f"COMBINED RESULTS ({n_total:,} total realisations)")
    print(f"{'='*60}")
    print(f"  R_NULL_95  (95th percentile):  {R_NULL_95:.6f}")
    print(f"  R_NULL_99  (99th percentile):  {R_NULL_99:.6f}")
    print(f"  R_NULL_995 (99.5th percentile): {R_NULL_995:.6f}")
    print(f"\n  *** R_NULL_99 = {R_NULL_99:.6f} ***")
    print(f"  This value goes into manifest.py before sealing.")
    
    results = {
        "version": "BOAT_v2.0_MC",
        "n_surveys_completed": sum(1 for v in per_survey_results.values() if v["status"] == "completed"),
        "n_total_realisations": n_total,
        "mc_base_seed": mc_base_seed,
        "decoy_seed": DECOY_AXES_SEED,
        "n_decoy_axes": N_DECOY_AXES,
        "r_weight_deg": R_WEIGHT_DEG,
        "combined_percentiles": {
            "R_NULL_95": round(R_NULL_95, 6),
            "R_NULL_99": round(R_NULL_99, 6),
            "R_NULL_995": round(R_NULL_995, 6),
        },
        "combined_stats": {
            "min": round(float(all_ratios.min()), 6),
            "max": round(float(all_ratios.max()), 6),
            "mean": round(float(all_ratios.mean()), 6),
            "std": round(float(all_ratios.std()), 6),
        },
        "per_survey": per_survey_results,
    }
    
    return results


def save_mc_results(results, filepath="boat/results/montecarlo_null_calibration.json"):
    """Save Monte Carlo results to JSON with SHA-256 hash."""
    from boat.hashing import hash_dict
    
    results["sha256"] = hash_dict(results)
    
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {filepath}")
    print(f"SHA-256: {results['sha256']}")
    
    return filepath


# ============================================================
# CONVENIENCE ENTRY POINT
# ============================================================

if __name__ == "__main__":
    # Default: run on whatever surveys are available
    # For testing, use n_realisations=100 first
    results = run_null_calibration(n_realisations=100)
    save_mc_results(results)
```

---

## IMPORTANT IMPLEMENTATION NOTES

### 1. Survey loaders are intentionally incomplete

The `load_6dfgs()` and `load_sdss()` functions contain file-inspection code and raise `NotImplementedError`. This is deliberate. **You need to:**

1. Run the file-inspection code to see the actual column names and delimiters
2. Complete the loaders based on what you find
3. Ensure each returns a DataFrame with exactly three columns: `'ra'`, `'dec'`, `'z'`

The 6dFGS file is typically pipe-delimited (`|`) with columns like `RAJ2000`, `DEJ2000`, `zfinal` or similar. The SDSS file is a CasJobs CSV export.

**Do NOT guess column names.** Inspect the actual files.

### 2. hash_dict may not exist in hashing.py

The v1.1 hashing.py has `hash_file` and `hash_dataframe`. It may not have `hash_dict`. If `hash_dict` doesn't exist, add it to hashing.py:

```python
def hash_dict(d):
    """SHA-256 hash of a dictionary (JSON-serialised)."""
    import hashlib
    import json
    # Remove any existing sha256 key to avoid circular hashing
    d_clean = {k: v for k, v in d.items() if k != 'sha256'}
    json_str = json.dumps(d_clean, sort_keys=True, default=str)
    return hashlib.sha256(json_str.encode()).hexdigest()
```

### 3. Testing sequence

After building the module:

**Test 1 — Smoke test (fast, 10 realisations):**
```python
from boat.montecarlo import run_null_calibration
results = run_null_calibration(surveys=["6dFGS"], n_realisations=10)
print(f"Smoke test R_NULL_99: {results['combined_percentiles']['R_NULL_99']}")
```
This should complete in under 2 minutes. It verifies the full pipeline works end-to-end.

**Test 2 — Check null distribution is centred around ~1.0:**
The ratio |r_true| / mean(|r_decoy|) under the null should have a mean near 1.0 (true axis is just another random direction when there's no signal). If the mean is significantly different from 1.0, something is wrong.

**Test 3 — Small calibration run (100 realisations per survey):**
```python
results = run_null_calibration(surveys=["6dFGS", "SDSS"], n_realisations=100)
```
This gives a rough R_NULL_99 estimate. The final production run will use 10,000.

---

## REPORT BACK

Once the module is built and Test 1 passes, report:

1. 6dFGS file format: what columns did you find? What delimiter?
2. SDSS file format: what columns did you find?
3. Survey loader status: both complete and returning correct DataFrame
4. hash_dict: already existed / added to hashing.py
5. Smoke test (10 realisations on 6dFGS): 
   - Completed successfully: yes/no
   - Time taken
   - Mean null ratio (should be ~1.0)
   - Any warnings or errors
6. Footprint weight computation: completed, time taken, weight statistics

**Do not run the full 10,000-realisation calibration yet.** We will do that as a separate step after verifying the pipeline.

Wait for the next instruction from Graham.

---

*Build prompt prepared by Claude Opus (Anthropic), 21 February 2026*
