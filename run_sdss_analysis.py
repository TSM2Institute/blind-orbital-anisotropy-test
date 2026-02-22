"""Run BOAT v2.0 test on 6dFGS + SDSS via runner.py."""

from boat.runner import run_production, save_production_results

results = run_production(survey_list=["6dFGS", "SDSS"], label="test_dev")
save_production_results(results)

print(f"\nOverall verdict: {results['overall_verdict']}")
print(f"Passes: {results['n_actual_passes']} of {results['n_qualified_surveys']} "
      f"(required: {results['n_required_passes']})")
