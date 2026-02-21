"""
Run 6dFGS full 10,000-realisation batch Monte Carlo calibration.
"""

from boat.montecarlo import run_mc_batch

result_6dfgs = run_mc_batch("6dFGS", batch_size=1000, total_realisations=10000)
print(f"\n6dFGS complete: {result_6dfgs['n_realisations']:,} realisations")
print(f"6dFGS p99: {result_6dfgs['percentiles']['p99']:.6f}")
print(f"6dFGS mean ratio: {result_6dfgs['mean_ratio']:.6f}")
print("\nDONE")
