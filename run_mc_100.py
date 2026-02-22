"""
Run SDSS full 10,000-realisation batch Monte Carlo calibration.
"""

from boat.montecarlo import run_mc_batch

result_sdss = run_mc_batch("SDSS", batch_size=1000, total_realisations=10000)
print(f"\nSDSS complete: {result_sdss['n_realisations']:,} realisations")
print(f"SDSS p99: {result_sdss['percentiles']['p99']:.6f}")
print(f"SDSS mean ratio: {result_sdss['mean_ratio']:.6f}")
print("\nDONE")
