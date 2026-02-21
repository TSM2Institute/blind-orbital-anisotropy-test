# BOAT v2.0 — Build Prompt #3b: Batch Monte Carlo Calibration

**Context:** montecarlo.py is built and verified. 100-realisation test passed on both 6dFGS and SDSS. We now run the full 10,000 realisations per survey, batched in 1,000-realisation chunks to survive Replit session limits.

---

## WHAT WE'RE DOING

Running 10,000 null realisations per survey (6dFGS and SDSS), saving progress after each 1,000-realisation batch. If a session dies, we resume from the last completed batch.

**Estimated time:**
- 6dFGS: ~80 seconds per 1,000 realisations × 10 batches = ~13 minutes total
- SDSS: ~3,700 seconds per 1,000 realisations × 10 batches = ~10.3 hours total

We run 6dFGS fully first (fast), then SDSS in batches.

---

## TASK: ADD BATCH FUNCTIONS TO montecarlo.py

Add the following functions to the END of `boat/montecarlo.py` (before the `if __name__` block). Do not modify any existing functions.

```python
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
    
    # Check for existing progress
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
    
    # Survey index for seed offset
    survey_names = list(SURVEY_LOADERS.keys())
    s_idx = survey_names.index(survey_name)
    
    # Load and prepare data (once)
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
    
    # Precompute (once)
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
    
    # Save metadata (for resume)
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
    
    # Run remaining batches
    all_ratios = list(existing_ratios)
    remaining = total_realisations - completed
    n_batches = int(np.ceil(remaining / batch_size))
    
    for b in range(n_batches):
        batch_start = completed + b * batch_size
        this_batch = min(batch_size, total_realisations - batch_start)
        
        print(f"\n  --- Batch {b+1}/{n_batches}: realisations {batch_start+1}–{batch_start+this_batch} ---")
        
        # Seed for this batch: deterministic from base seed + survey + batch offset
        batch_seed = mc_base_seed + s_idx * 100000 + batch_start
        
        t_batch = time.time()
        ratios = run_mc_for_survey(
            ra, dec, z_obs, weights,
            cos_theta_true, cos_theta_decoys,
            this_batch, batch_seed
        )
        elapsed = time.time() - t_batch
        
        all_ratios.extend(ratios.tolist())
        
        # Save progress after each batch
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
        # Auto-detect completed surveys
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
```

**IMPORTANT implementation notes:**

1. The batch seeding is deterministic: `batch_seed = mc_base_seed + s_idx * 100000 + batch_start`. This means resuming produces identical results to a continuous run — each batch gets a unique, reproducible seed.

2. Progress is saved as `.npy` (binary numpy) for the ratio arrays and `.json` for metadata. This is in `boat/results/mc_progress/` — separate from the final production results.

3. `combine_batch_results()` is called AFTER all surveys complete to produce the final combined output.

---

## EXECUTION PLAN

### Step 1: Run 6dFGS full calibration (fast — ~13 minutes)

```python
from boat.montecarlo import run_mc_batch

result_6dfgs = run_mc_batch("6dFGS", batch_size=1000, total_realisations=10000)
print(f"\n6dFGS complete: {result_6dfgs['n_realisations']:,} realisations")
print(f"6dFGS p99: {result_6dfgs['percentiles']['p99']:.6f}")
```

Report the result, then wait for instruction before starting SDSS.

### Step 2 (separate session if needed): Run SDSS in batches

```python
from boat.montecarlo import run_mc_batch

result_sdss = run_mc_batch("SDSS", batch_size=1000, total_realisations=10000)
print(f"\nSDSS complete: {result_sdss['n_realisations']:,} realisations")
print(f"SDSS p99: {result_sdss['percentiles']['p99']:.6f}")
```

If the session dies mid-run, just paste the same command again — it will resume from the last saved batch.

### Step 3: Combine and save final results

```python
from boat.montecarlo import combine_batch_results, save_mc_results

results = combine_batch_results()
save_mc_results(results, filepath="boat/results/montecarlo_null_calibration_partial.json")
```

(Named "partial" because DESI surveys will be added later.)

---

## REPORT BACK AFTER STEP 1 (6dFGS only)

1. All 10 batches completed: yes/no
2. Total time for 6dFGS
3. 6dFGS null ratio percentiles (p50, p95, p99, p99.5)
4. 6dFGS mean null ratio (should be ~1.0)
5. Progress files saved in boat/results/mc_progress/

**Then wait for instruction before starting SDSS.**

---

*Build prompt prepared by Claude Opus (Anthropic), 21 February 2026*
