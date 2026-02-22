"""DESI BGS MC calibration — 1,000 realisations on 1M subsample."""

import sys
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

from boat.datasets import load_survey, apply_quality_cuts, check_prequalification
from boat.montecarlo import subsample_survey, run_mc_batch

print("=" * 60)
print("DESI BGS — Subsample Validation Run")
print("=" * 60)

df = load_survey("DESI_BGS")
df = apply_quality_cuts(df, "DESI_BGS")
passes, sky_frac, med_z = check_prequalification(df, "DESI_BGS")
print(f"Full BGS: {len(df):,} galaxies, sky={sky_frac:.1%}, median_z={med_z:.4f}")

df_sub = subsample_survey(df, max_galaxies=1_000_000, seed=12345)

print(f"\nSubsample footprint check:")
ra_bins_full = np.floor(df['ra'].values).astype(int)
dec_bins_full = np.floor(df['dec'].values + 90).astype(int)
bins_full = len(np.unique(ra_bins_full * 180 + dec_bins_full))

ra_bins_sub = np.floor(df_sub['ra'].values).astype(int)
dec_bins_sub = np.floor(df_sub['dec'].values + 90).astype(int)
bins_sub = len(np.unique(ra_bins_sub * 180 + dec_bins_sub))

print(f"  Full: {bins_full} occupied sky bins")
print(f"  Subsample: {bins_sub} occupied sky bins")
print(f"  Retention: {bins_sub/bins_full:.1%} of footprint preserved")

df_sub.to_csv("attached_assets/DESI_DR1_BGS_MC_SUBSAMPLE.csv", index=False)
print(f"\nSubsample saved. Running 1,000 MC realisations...")

t0 = time.time()
result = run_mc_batch("DESI_BGS_MC", batch_size=1000, total_realisations=1000)
elapsed = time.time() - t0

print(f"\n{'=' * 60}")
print(f"DESI BGS (1M subsample) — 1,000 realisations:")
print(f"  Total time: {elapsed:.1f}s")
print(f"  Mean ratio: {result['mean_ratio']:.4f}")
print(f"  p95: {result['percentiles']['p95']:.6f}")
print(f"  p99: {result['percentiles']['p99']:.6f}")
print(f"  p99.5: {result['percentiles']['p995']:.6f}")

print(f"\nNull distribution comparison (p99 values):")
print(f"  6dFGS (120k, full):        p99 = 1.5776")
print(f"  SDSS (1.58M, full):        p99 = 1.4984")
print(f"  DESI BGS (1M subsample):   p99 = {result['percentiles']['p99']:.4f}")
print(f"\nIf BGS p99 is within ~0.1 of the others, subsampling is validated.")
print("=" * 60)
