"""DESI MC calibration — LRG subsample creation + full 10,000 runs for BGS and LRG."""

import sys
import time

sys.stdout.reconfigure(line_buffering=True)

from boat.datasets import load_survey, apply_quality_cuts
from boat.montecarlo import subsample_survey, run_mc_batch, combine_batch_results, save_mc_results

print("=" * 60)
print("STEP 1: Create LRG subsample")
print("=" * 60)

df = load_survey("DESI_LRG")
df = apply_quality_cuts(df, "DESI_LRG")
print(f"Full LRG after cuts: {len(df):,} galaxies")

df_sub = subsample_survey(df, max_galaxies=1_000_000, seed=12345)
df_sub.to_csv("attached_assets/DESI_DR1_LRG_MC_SUBSAMPLE.csv", index=False)
print(f"LRG subsample saved: {len(df_sub):,} galaxies")
del df, df_sub

print(f"\n{'=' * 60}")
print("STEP 2: Run full 10,000 MC on BGS subsample")
print("=" * 60)

t0_bgs = time.time()
result_bgs = run_mc_batch("DESI_BGS_MC", batch_size=1000, total_realisations=10000)
elapsed_bgs = time.time() - t0_bgs

print(f"\nDESI BGS MC complete:")
print(f"  Total time: {elapsed_bgs:.1f}s ({elapsed_bgs/60:.1f} min)")
print(f"  Mean ratio: {result_bgs['mean_ratio']:.4f}")
print(f"  p99: {result_bgs['percentiles']['p99']:.6f}")

print(f"\n{'=' * 60}")
print("STEP 3: Run full 10,000 MC on LRG subsample")
print("=" * 60)

t0_lrg = time.time()
result_lrg = run_mc_batch("DESI_LRG_MC", batch_size=1000, total_realisations=10000)
elapsed_lrg = time.time() - t0_lrg

print(f"\nDESI LRG MC complete:")
print(f"  Total time: {elapsed_lrg:.1f}s ({elapsed_lrg/60:.1f} min)")
print(f"  Mean ratio: {result_lrg['mean_ratio']:.4f}")
print(f"  p99: {result_lrg['percentiles']['p99']:.6f}")

print(f"\n{'=' * 60}")
print("STEP 4: Combine all 4 surveys → final R_NULL_99")
print("=" * 60)

results = combine_batch_results(
    survey_names=["6dFGS", "SDSS", "DESI_BGS_MC", "DESI_LRG_MC"]
)
save_mc_results(results, filepath="boat/results/montecarlo_null_calibration_final.json")

print(f"\n{'=' * 60}")
print(f"FINAL COMBINED R_NULL_99 = {results['combined_percentiles']['R_NULL_99']:.6f}")
print(f"From {results['n_total_realisations']:,} total realisations across 4 surveys")
print(f"{'=' * 60}")

print(f"\nPer-survey p99 comparison:")
print(f"  {'Survey':<25} {'N galaxies':>12} {'Realisations':>14} {'p99':>10}")
print(f"  {'-'*25} {'-'*12} {'-'*14} {'-'*10}")
for name, s in results['per_survey'].items():
    print(f"  {name:<25} {s['n_galaxies']:>12,} {s['n_realisations']:>14,} {s['percentiles']['p99']:>10.4f}")

print(f"\n  {'Combined':.<25} {'':>12} {results['n_total_realisations']:>14,} {results['combined_percentiles']['R_NULL_99']:>10.4f}")
print(f"\nTotal compute time: {(elapsed_bgs + elapsed_lrg):.0f}s ({(elapsed_bgs + elapsed_lrg)/3600:.1f} hours)")
