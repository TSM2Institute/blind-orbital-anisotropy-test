"""
BOAT v2.0 — Production Runner

Executes the complete blind test protocol:
1. Load and pre-qualify all available surveys
2. Run full analysis on each qualifying survey
3. Evaluate per-dataset PASS/FAIL
4. Compute overall verdict
5. Save all results with SHA-256 hashes

Usage:
    python -m boat.runner          # production run
    python -m boat.runner --test   # test run (6dFGS only)
"""

import numpy as np
import json
import math
import time
import sys
import os

from boat.manifest import (
    R_NULL_99, N_DECOY_AXES, DECOY_AXES_SEED, R_WEIGHT_DEG,
    PASS_REQUIRED_FORMULA, P_VALUE_THRESHOLD, RANK_THRESHOLD,
    V_HAT, R_HAT, K_HAT
)
from boat.datasets import (
    load_survey, apply_quality_cuts, check_prequalification,
    get_available_surveys
)
from boat.analysis import run_full_analysis
from boat.hashing import hash_dataframe, hash_dict


def compute_pass_required(n_surveys):
    """
    Compute number of passes required: ceil(N/2) + 1

    N=3 → 3 (unanimous)
    N=4 → 3
    N=5 → 4
    """
    return math.ceil(n_surveys / 2) + 1


def run_production(survey_list=None, label="production"):
    """
    Run the full BOAT v2.0 blind test.

    Args:
        survey_list: list of survey names (default: all available)
        label: run label for output filename

    Returns:
        dict with complete results
    """
    t_start = time.time()

    print("=" * 70)
    print("BOAT v2.0 — BLIND ORBITAL ANISOTROPY TEST")
    print(f"R_NULL_99 = {R_NULL_99}")
    print("=" * 70)

    if R_NULL_99 is None:
        print("\nWARNING: R_NULL_99 is not set. Running in test mode.")
        print("All ratio tests will FAIL. This is expected during development.")

    if survey_list is None:
        survey_list = get_available_surveys()

    print(f"\nSurveys to process: {survey_list}")

    # ============================================================
    # PHASE 1: Load, cut, and pre-qualify each survey
    # ============================================================

    qualified_surveys = []
    disqualified_surveys = []
    survey_data = {}

    for survey_name in survey_list:
        print(f"\n{'='*60}")
        print(f"Loading: {survey_name}")
        print(f"{'='*60}")

        try:
            df = load_survey(survey_name)
        except (FileNotFoundError, NotImplementedError) as e:
            print(f"  SKIPPED: {e}")
            disqualified_surveys.append({
                "survey": survey_name,
                "reason": str(e),
            })
            continue

        print(f"  Raw count: {len(df):,}")

        df = apply_quality_cuts(df, survey_name)

        passes, sky_frac, med_z = check_prequalification(df, survey_name)

        if not passes:
            disqualified_surveys.append({
                "survey": survey_name,
                "reason": "Failed pre-qualification",
                "sky_fraction": sky_frac,
                "median_z": med_z,
                "n_galaxies": len(df),
            })
            continue

        data_hash = hash_dataframe(df)

        qualified_surveys.append(survey_name)
        survey_data[survey_name] = {
            "df": df,
            "sky_fraction": sky_frac,
            "median_z": med_z,
            "data_hash": data_hash,
        }

    n_qualified = len(qualified_surveys)
    print(f"\n{'='*60}")
    print(f"Qualified surveys: {n_qualified} — {qualified_surveys}")
    if disqualified_surveys:
        print(f"Disqualified: {[d['survey'] for d in disqualified_surveys]}")
    print(f"{'='*60}")

    if n_qualified == 0:
        print("\nERROR: No surveys qualified. Cannot proceed.")
        return {"error": "No surveys qualified"}

    n_required = compute_pass_required(n_qualified)
    print(f"Required passes for overall PASS: {n_required} of {n_qualified}")

    # ============================================================
    # PHASE 2: Run analysis on each qualified survey
    # ============================================================

    per_survey_results = {}
    n_passed = 0

    for survey_name in qualified_surveys:
        sd = survey_data[survey_name]
        df = sd["df"]

        results = run_full_analysis(
            df['ra'].values,
            df['dec'].values,
            df['z'].values,
            survey_name
        )

        results["sky_fraction"] = sd["sky_fraction"]
        results["median_z"] = sd["median_z"]
        results["data_hash"] = sd["data_hash"]

        per_survey_results[survey_name] = results

        if results["verdict"] == "PASS":
            n_passed += 1

    # ============================================================
    # PHASE 3: Overall verdict
    # ============================================================

    overall_pass = n_passed >= n_required

    orbital_detections = []
    for sname, res in per_survey_results.items():
        od = res.get("orbital_discriminator", {})
        if od.get("detection", False):
            orbital_detections.append(sname)

    if orbital_detections:
        orbital_summary = f"DETECTED in {len(orbital_detections)} of {n_qualified} surveys: {orbital_detections}"
    else:
        orbital_summary = "NOT DETECTED in any survey"

    elapsed = time.time() - t_start

    # ============================================================
    # COMPILE FINAL RESULTS
    # ============================================================

    print(f"\n{'='*70}")
    print(f"BOAT v2.0 — FINAL RESULTS")
    print(f"{'='*70}")
    print(f"  Qualified surveys: {n_qualified}")
    print(f"  Required passes: {n_required}")
    print(f"  Actual passes: {n_passed}")
    print(f"  Overall verdict: {'PASS' if overall_pass else 'FAIL'}")
    print(f"  Orbital discriminator: {orbital_summary}")
    print(f"  Total runtime: {elapsed:.1f}s")

    for sname, res in per_survey_results.items():
        pt = res["primary_test"]
        print(f"\n  {sname}:")
        print(f"    Verdict: {res['verdict']}")
        print(f"    |r| = {pt['true_axis']['abs_r']:.6f}, "
              f"ratio = {pt['ratio']:.4f}, "
              f"rank = #{pt['rank']}, "
              f"p = {pt['true_axis']['p_value']:.2e}")
        print(f"    z-bin: {'PASS' if res['zbin_invariance']['pass'] else 'FAIL'}")
        print(f"    Orbital: {res['orbital_discriminator']['verdict']}")

    print(f"\n{'='*70}")

    final_results = {
        "version": "BOAT_v2.0",
        "label": label,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "R_NULL_99": R_NULL_99,
        "n_decoy_axes": N_DECOY_AXES,
        "decoy_seed": DECOY_AXES_SEED,
        "r_weight_deg": R_WEIGHT_DEG,
        "n_qualified_surveys": n_qualified,
        "n_required_passes": n_required,
        "n_actual_passes": n_passed,
        "overall_verdict": "PASS" if overall_pass else "FAIL",
        "orbital_discriminator_summary": orbital_summary,
        "elapsed_seconds": round(elapsed, 1),
        "qualified_surveys": qualified_surveys,
        "disqualified_surveys": disqualified_surveys,
        "per_survey": per_survey_results,
    }

    return final_results


def save_production_results(results, filepath=None):
    """Save results JSON with SHA-256 hash."""
    if filepath is None:
        label = results.get("label", "production")
        filepath = f"boat/results/boat_results_{label}_v2.0.json"

    results["sha256"] = hash_dict(results)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.bool_):
                return bool(obj)
            return super().default(obj)

    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)

    print(f"\nResults saved to: {filepath}")
    print(f"SHA-256: {results['sha256']}")

    return filepath


# ============================================================
# COMMAND LINE INTERFACE
# ============================================================

if __name__ == "__main__":
    if "--test" in sys.argv:
        print("TEST MODE: Running on 6dFGS only")
        results = run_production(survey_list=["6dFGS"], label="test")
    else:
        print("PRODUCTION MODE: Running on all available surveys")
        results = run_production(label="production")

    save_production_results(results)
