"""
Run 100-realisation Monte Carlo calibration on both surveys.
Saves results to boat/results/mc_calibration_100.json
"""

from boat.montecarlo import run_null_calibration, save_mc_results

results = run_null_calibration(
    surveys=["6dFGS", "SDSS"],
    n_realisations=100
)

save_mc_results(results, filepath="boat/results/mc_calibration_100.json")

print("\n\nDONE — results saved to boat/results/mc_calibration_100.json")
