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
import pandas as pd
import json
import time
import os

from boat.manifest import (
    V_OBS, C_KM_S, V_HAT, R_HAT, K_HAT, V_OBS_VEC,
    N_DECOY_AXES, DECOY_AXES_SEED, R_WEIGHT_DEG,
    N_MC_REALISATIONS
)
from boat.datasets import (
    load_6dfgs, load_sdss, apply_quality_cuts,
    check_prequalification, SURVEY_LOADERS, LRG_SURVEYS
)


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
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    unit_vecs = np.column_stack([x, y, z])

    chord_threshold = 2.0 * np.sin(np.radians(R_WEIGHT_DEG) / 2.0)

    tree = cKDTree(unit_vecs)
    counts = tree.query_ball_point(unit_vecs, r=chord_threshold, return_length=True)

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
        z_shuffled = rng.permutation(z_obs)

        dz = z_shuffled - np.mean(z_shuffled)

        r_true = weighted_pearson_r(cos_theta_true, dz, weights)

        r_decoys = np.empty(n_decoys)
        for j in range(n_decoys):
            r_decoys[j] = weighted_pearson_r(cos_theta_decoys[j], dz, weights)

        mean_decoy_abs_r = np.mean(np.abs(r_decoys))
        if mean_decoy_abs_r == 0:
            ratios[i] = 0.0
        else:
            ratios[i] = np.abs(r_true) / mean_decoy_abs_r

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

    decoy_axes = generate_decoy_axes()
    true_axis = V_HAT.reshape(1, 3)
    all_axes = np.vstack([true_axis, decoy_axes])

    print(f"\nDecoy axes generated: {len(decoy_axes)} (seed: {DECOY_AXES_SEED})")
    print(f"Total axes per dataset: {len(all_axes)}")

    all_ratios = []
    per_survey_results = {}

    for s_idx, survey_name in enumerate(surveys):
        print(f"\n{'='*60}")
        print(f"Processing: {survey_name}")
        print(f"{'='*60}")

        t_start = time.time()

        loader = SURVEY_LOADERS[survey_name]
        df = loader()
        print(f"  Loaded: {len(df):,} galaxies")

        df = apply_quality_cuts(df, survey_name)

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

        print(f"\n  Computing footprint weights (R = {R_WEIGHT_DEG}°)...")
        t_weights = time.time()
        weights = compute_footprint_weights(ra, dec)
        print(f"  Weights computed in {time.time() - t_weights:.1f}s")
        print(f"  Weight stats: min={weights.min():.4f}, max={weights.max():.4f}, "
              f"mean={weights.mean():.4f}, std={weights.std():.4f}")

        unit_vecs = ra_dec_to_unit_vectors(ra, dec)
        cos_theta_all = precompute_cos_theta(unit_vecs, all_axes)
        cos_theta_true = cos_theta_all[0]
        cos_theta_decoys = cos_theta_all[1:]

        print(f"\n  Running {n_realisations:,} null realisations...")
        survey_seed = mc_base_seed + s_idx
        ratios = run_mc_for_survey(
            ra, dec, z_obs, weights,
            cos_theta_true, cos_theta_decoys,
            n_realisations, survey_seed
        )

        elapsed = time.time() - t_start

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
# BATCH MONTE CARLO RUNNER
# ============================================================

def run_mc_batch(survey_name, batch_size=1000, total_realisations=10000, mc_base_seed=42):
    """
    Run Monte Carlo in batches, saving progress after each batch.

    Precomputes weights and cos(θ) once, then runs batch_size
    realisations at a time, appending ratios to a file.

    Resume-safe: checks for existing partial results and continues
    from where it left off.

    Args:
        survey_name: "6dFGS" or "SDSS"
        batch_size: realisations per batch (default 1000)
        total_realisations: target total (default 10000)
        mc_base_seed: base seed for RNG

    Returns:
        dict with survey results
    """
    import json

    progress_dir = "boat/results/mc_progress"
    os.makedirs(progress_dir, exist_ok=True)
    ratios_file = os.path.join(progress_dir, f"{survey_name}_ratios.npy")
    meta_file = os.path.join(progress_dir, f"{survey_name}_meta.json")

    existing_ratios = np.array([])
    if os.path.exists(ratios_file) and os.path.exists(meta_file):
        existing_ratios = np.load(ratios_file)
        with open(meta_file, 'r') as f:
            meta = json.load(f)
        completed = len(existing_ratios)
        print(f"  Resuming {survey_name}: {completed:,} realisations already done")
        if completed >= total_realisations:
            print(f"  {survey_name} already complete ({completed:,}/{total_realisations:,})")
            return _compile_survey_results(survey_name, existing_ratios, meta)
    else:
        completed = 0

    survey_names = list(SURVEY_LOADERS.keys())
    s_idx = survey_names.index(survey_name)

    print(f"\n{'='*60}")
    print(f"BATCH MC: {survey_name}")
    print(f"Target: {total_realisations:,} realisations in batches of {batch_size:,}")
    print(f"{'='*60}")

    loader = SURVEY_LOADERS[survey_name]
    df = loader()
    print(f"  Loaded: {len(df):,} galaxies")

    df = apply_quality_cuts(df, survey_name)
    passes, sky_frac, med_z = check_prequalification(df, survey_name)
    if not passes:
        print(f"  {survey_name} does not pre-qualify. Aborting.")
        return {"status": "disqualified", "survey": survey_name}

    ra = df['ra'].values
    dec = df['dec'].values
    z_obs = df['z'].values

    print(f"\n  Computing footprint weights...")
    t0 = time.time()
    weights = compute_footprint_weights(ra, dec)
    print(f"  Weights done in {time.time()-t0:.1f}s")

    decoy_axes = generate_decoy_axes()
    true_axis = V_HAT.reshape(1, 3)
    all_axes = np.vstack([true_axis, decoy_axes])

    unit_vecs = ra_dec_to_unit_vectors(ra, dec)
    cos_theta_all = precompute_cos_theta(unit_vecs, all_axes)
    cos_theta_true = cos_theta_all[0]
    cos_theta_decoys = cos_theta_all[1:]

    meta = {
        "survey": survey_name,
        "n_galaxies": len(df),
        "sky_fraction": float(sky_frac),
        "median_z": float(med_z),
        "mc_base_seed": mc_base_seed,
        "batch_size": batch_size,
        "total_target": total_realisations,
        "weight_stats": {
            "min": float(weights.min()),
            "max": float(weights.max()),
            "mean": float(weights.mean()),
            "std": float(weights.std()),
        }
    }

    all_ratios = list(existing_ratios)
    remaining = total_realisations - completed
    n_batches = int(np.ceil(remaining / batch_size))

    for b in range(n_batches):
        batch_start = completed + b * batch_size
        this_batch = min(batch_size, total_realisations - batch_start)

        print(f"\n  --- Batch {b+1}/{n_batches}: realisations {batch_start+1}–{batch_start+this_batch} ---")

        batch_seed = mc_base_seed + s_idx * 100000 + batch_start

        t_batch = time.time()
        ratios = run_mc_for_survey(
            ra, dec, z_obs, weights,
            cos_theta_true, cos_theta_decoys,
            this_batch, batch_seed
        )
        elapsed = time.time() - t_batch

        all_ratios.extend(ratios.tolist())

        np.save(ratios_file, np.array(all_ratios))
        with open(meta_file, 'w') as f:
            json.dump(meta, f, indent=2)

        done = batch_start + this_batch
        print(f"  Batch complete: {elapsed:.1f}s | "
              f"Progress: {done:,}/{total_realisations:,} | "
              f"Batch mean ratio: {np.mean(ratios):.4f}")

    all_ratios = np.array(all_ratios)
    return _compile_survey_results(survey_name, all_ratios, meta)


def _compile_survey_results(survey_name, ratios, meta):
    """Compile final per-survey results dict."""
    return {
        "status": "completed",
        "survey": survey_name,
        "n_galaxies": meta["n_galaxies"],
        "sky_fraction": meta["sky_fraction"],
        "median_z": meta["median_z"],
        "n_realisations": len(ratios),
        "mc_base_seed": meta["mc_base_seed"],
        "weight_stats": meta["weight_stats"],
        "percentiles": {
            "p50": round(float(np.percentile(ratios, 50)), 6),
            "p95": round(float(np.percentile(ratios, 95)), 6),
            "p99": round(float(np.percentile(ratios, 99)), 6),
            "p995": round(float(np.percentile(ratios, 99.5)), 6),
        },
        "min_ratio": round(float(ratios.min()), 6),
        "max_ratio": round(float(ratios.max()), 6),
        "mean_ratio": round(float(ratios.mean()), 6),
    }


def combine_batch_results(survey_names=None):
    """
    Combine completed batch results from all surveys into final
    calibration output.

    Call this after all surveys have completed their 10,000 runs.

    Returns:
        dict with combined R_NULL_99 and per-survey breakdowns
    """
    import json

    progress_dir = "boat/results/mc_progress"

    if survey_names is None:
        survey_names = []
        for f in os.listdir(progress_dir):
            if f.endswith("_ratios.npy"):
                name = f.replace("_ratios.npy", "")
                survey_names.append(name)

    print(f"Combining results from: {survey_names}")

    all_ratios = []
    per_survey = {}

    for name in survey_names:
        ratios_file = os.path.join(progress_dir, f"{name}_ratios.npy")
        meta_file = os.path.join(progress_dir, f"{name}_meta.json")

        if not os.path.exists(ratios_file):
            print(f"  WARNING: {name} has no results file")
            continue

        ratios = np.load(ratios_file)
        with open(meta_file, 'r') as f:
            meta = json.load(f)

        per_survey[name] = _compile_survey_results(name, ratios, meta)
        all_ratios.extend(ratios.tolist())

        print(f"  {name}: {len(ratios):,} realisations, "
              f"p99 = {per_survey[name]['percentiles']['p99']:.6f}")

    all_ratios = np.array(all_ratios)
    n_total = len(all_ratios)

    R_NULL_95 = float(np.percentile(all_ratios, 95))
    R_NULL_99 = float(np.percentile(all_ratios, 99))
    R_NULL_995 = float(np.percentile(all_ratios, 99.5))

    print(f"\n{'='*60}")
    print(f"COMBINED: {n_total:,} total realisations from {len(per_survey)} surveys")
    print(f"  R_NULL_95:  {R_NULL_95:.6f}")
    print(f"  R_NULL_99:  {R_NULL_99:.6f}")
    print(f"  R_NULL_995: {R_NULL_995:.6f}")
    print(f"\n  *** R_NULL_99 = {R_NULL_99:.6f} ***")
    print(f"{'='*60}")

    results = {
        "version": "BOAT_v2.0_MC",
        "n_surveys_completed": len(per_survey),
        "surveys_included": list(per_survey.keys()),
        "n_total_realisations": n_total,
        "mc_base_seed": 42,
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
        "per_survey": per_survey,
        "note": "Partial calibration (6dFGS + SDSS only). DESI surveys to be added after data download. Final R_NULL_99 will be recomputed from all 4 surveys before sealing."
    }

    return results


# ============================================================
# CONVENIENCE ENTRY POINT
# ============================================================

if __name__ == "__main__":
    results = run_null_calibration(n_realisations=100)
    save_mc_results(results)
