# BOAT v2.0 — Build Prompt #5: datasets.py + runner.py

**Context:** manifest.py, montecarlo.py, and analysis.py are all built and verified. 6dFGS and SDSS diagnostics complete. We now build the final two modules: the survey loading layer and the production runner.

---

## MODULE 1: boat/datasets.py

This module handles survey loading, quality cuts, and pre-qualification. It consolidates the loaders currently inside montecarlo.py into the canonical location, and adds stubs for DESI.

**Important:** montecarlo.py already has working `load_6dfgs()`, `load_sdss()`, `apply_quality_cuts()`, and `check_prequalification()`. You have two options:

**Option A (recommended):** Move the working implementations from montecarlo.py into datasets.py, then update montecarlo.py to import from datasets.py. This eliminates duplication.

**Option B:** Copy the implementations into datasets.py and leave montecarlo.py as-is. Slight duplication but zero risk of breaking the MC module.

**Choose Option A** — cleaner codebase, and montecarlo.py's functionality is verified so refactoring the import is low-risk.

### datasets.py structure:

```python
"""
BOAT v2.0 — Survey Data Loaders and Pre-Qualification

Loads, cleans, and pre-qualifies galaxy survey data for the
blind orbital anisotropy test.

Supported surveys:
- 6dFGS DR3 (Southern, available)
- SDSS DR18 (Northern, available)
- DESI DR1 BGS (Both hemispheres, pending data download)
- DESI DR1 LRG (Both hemispheres, pending data download)

Usage:
    from boat.datasets import load_survey, get_available_surveys
    df = load_survey("6dFGS")
"""

import numpy as np
import pandas as pd
import os

from boat.manifest import (
    Z_MIN, Z_MAX, Z_MIN_LRG, Z_MAX_LRG,
    GALACTIC_PLANE_MASK_DEG, MIN_SAMPLE_SIZE,
    MIN_SKY_FRACTION, MIN_MEDIAN_Z,
    DUPLICATE_MATCH_ARCSEC, DATA_DIR, DATA_FILES,
    C_KM_S
)

# Which surveys use LRG redshift range
LRG_SURVEYS = {"DESI_LRG"}


def load_6dfgs():
    """
    Load 6dFGS DR3 survey data.
    
    [MOVE the working implementation from montecarlo.py here]
    
    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    # ... (move from montecarlo.py)


def load_sdss():
    """
    Load SDSS DR18 survey data.
    
    [MOVE the working implementation from montecarlo.py here]
    
    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    # ... (move from montecarlo.py)


def load_desi_bgs():
    """
    Load DESI DR1 BGS survey data.
    
    Data file not yet available — will be downloaded from
    NOIRLab Astro Data Lab after account approval.
    
    Expected file: attached_assets/DESI_DR1_BGS.csv (or .fits)
    Expected columns: RA, DEC, Z (spectroscopic)
    Quality flag: ZWARN == 0
    
    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, "DESI_DR1_BGS.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"DESI BGS data not found at {filepath}. "
            "Download from NOIRLab Astro Data Lab first."
        )
    # TODO: Complete when data format is known
    raise NotImplementedError("DESI BGS loader — complete after data download")


def load_desi_lrg():
    """
    Load DESI DR1 LRG survey data.
    
    Same source as BGS but different target class and redshift range.
    
    Expected file: attached_assets/DESI_DR1_LRG.csv (or .fits)
    
    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, "DESI_DR1_LRG.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"DESI LRG data not found at {filepath}. "
            "Download from NOIRLab Astro Data Lab first."
        )
    # TODO: Complete when data format is known
    raise NotImplementedError("DESI LRG loader — complete after data download")


# Registry of available loaders
SURVEY_LOADERS = {
    "6dFGS": load_6dfgs,
    "SDSS": load_sdss,
    "DESI_BGS": load_desi_bgs,
    "DESI_LRG": load_desi_lrg,
}


def get_available_surveys():
    """
    Return list of surveys whose data files are present.
    Does NOT check pre-qualification — just file existence.
    """
    available = []
    # Map survey names to expected files
    file_checks = {
        "6dFGS": os.path.join(DATA_DIR, DATA_FILES.get("6dFGS", "")),
        "SDSS": os.path.join(DATA_DIR, DATA_FILES.get("SDSS", "")),
        "DESI_BGS": os.path.join(DATA_DIR, "DESI_DR1_BGS.csv"),
        "DESI_LRG": os.path.join(DATA_DIR, "DESI_DR1_LRG.csv"),
    }
    for name, path in file_checks.items():
        if path and os.path.exists(path):
            available.append(name)
    return available


def load_survey(survey_name):
    """
    Load a survey by name.
    
    Args:
        survey_name: one of "6dFGS", "SDSS", "DESI_BGS", "DESI_LRG"
    
    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    if survey_name not in SURVEY_LOADERS:
        raise ValueError(f"Unknown survey: {survey_name}. "
                        f"Available: {list(SURVEY_LOADERS.keys())}")
    return SURVEY_LOADERS[survey_name]()


def apply_quality_cuts(df, survey_name):
    """
    Apply standard quality cuts to survey data.
    
    [MOVE the working implementation from montecarlo.py here]
    
    1. Redshift range (standard or LRG)
    2. Galactic plane mask (|b| >= 10°)
    3. Duplicate removal (match_coordinates_sky, < 1 arcsec)
    4. Minimum sample size check
    
    Returns:
        DataFrame after cuts
    """
    # ... (move from montecarlo.py)


def check_prequalification(df, survey_name):
    """
    Check if survey meets pre-qualification criteria.
    
    [MOVE the working implementation from montecarlo.py here]
    
    - Sky fraction >= 15%
    - Median z >= 0.04
    - N >= 10,000
    
    Returns:
        (passes: bool, sky_fraction: float, median_z: float)
    """
    # ... (move from montecarlo.py)
```

### After moving functions to datasets.py, update montecarlo.py imports:

Replace the local loader functions and quality cut functions in montecarlo.py with imports:

```python
from boat.datasets import (
    load_6dfgs, load_sdss, apply_quality_cuts, 
    check_prequalification, SURVEY_LOADERS, LRG_SURVEYS
)
```

Remove the duplicated function bodies from montecarlo.py but keep the import references so it works identically.

### Verify after refactoring:

```python
# Quick check that montecarlo still works after refactor
from boat.montecarlo import run_null_calibration
results = run_null_calibration(surveys=["6dFGS"], n_realisations=5)
print(f"MC smoke test: mean ratio = {results['combined_stats']['mean']:.4f}")
# Should complete without errors
```

---

## MODULE 2: boat/runner.py

The production orchestrator. Runs the complete BOAT v2.0 protocol on all qualifying surveys and produces the final results JSON.

```python
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
    
    # Determine which surveys to run
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
        
        # Quality cuts
        df = apply_quality_cuts(df, survey_name)
        
        # Pre-qualification
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
        
        # Hash the filtered data for reproducibility
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
        
        # Add pre-qualification info
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
    
    # Orbital discriminator summary
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
    
    # Remove DataFrames and large arrays before serialisation
    # (they're already summarised in the results)
    
    results["sha256"] = hash_dict(results)
    
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Custom serialiser for numpy types
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
```

---

## IMPLEMENTATION NOTES

### 1. hash_dict import

runner.py imports `hash_dict` from hashing.py. If this was already added during the montecarlo.py build, it will work. If not, verify it exists.

### 2. NumpyEncoder

The JSON serialiser needs to handle numpy types (int64, float64, bool_) that come from the analysis results. The custom encoder handles this.

### 3. The runner does NOT modify any manifest values

R_NULL_99 is read from manifest.py at import time. The runner never writes to manifest.py. The only way R_NULL_99 gets set is by manually editing manifest.py after Monte Carlo calibration.

### 4. DESI surveys will fail gracefully

When DESI data files aren't present, `load_survey("DESI_BGS")` raises `FileNotFoundError`, which the runner catches and reports as "SKIPPED". The overall verdict adjusts to the number of qualifying surveys.

---

## VERIFICATION

After building both modules, run this test:

```python
from boat.runner import run_production, save_production_results

# Test run on available surveys (6dFGS + SDSS)
results = run_production(survey_list=["6dFGS", "SDSS"], label="test_dev")
save_production_results(results)

print(f"\nOverall verdict: {results['overall_verdict']}")
print(f"Passes: {results['n_actual_passes']} of {results['n_qualified_surveys']} "
      f"(required: {results['n_required_passes']})")
```

**Expected outcome:**
- Both surveys qualify (sky >= 15%, median z >= 0.04)
- Both surveys FAIL (R_NULL_99 is None → ratio test fails automatically)
- Overall: FAIL
- All numbers populated, results JSON saved with SHA-256

Also verify montecarlo.py still works after the refactor:

```python
from boat.montecarlo import run_null_calibration
results = run_null_calibration(surveys=["6dFGS"], n_realisations=5)
print(f"MC smoke test: {results['combined_stats']['mean']:.4f}")
```

---

## REPORT BACK

1. datasets.py created: yes/no
2. montecarlo.py refactored (imports from datasets.py): yes/no
3. runner.py created: yes/no
4. MC smoke test after refactor: PASS/FAIL
5. Test run (6dFGS + SDSS):
   - Both surveys qualified: yes/no
   - Overall verdict: FAIL (expected)
   - Results JSON saved: yes/no
   - Any errors

**After this, the core pipeline is complete.** We wait for DESI data to finish the Monte Carlo calibration and run production.

---

*Build prompt prepared by Claude Opus (Anthropic), 22 February 2026*
