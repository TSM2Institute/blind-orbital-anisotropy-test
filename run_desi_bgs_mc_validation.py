"""BOAT v2.0 — PRODUCTION RUN. Results published regardless of outcome."""

import sys

sys.stdout.reconfigure(line_buffering=True)

from boat.runner import run_production, save_production_results

results = run_production(label="production")
filepath = save_production_results(results)

print(f"\n{'='*70}")
print(f"PRODUCTION RUN COMPLETE")
print(f"Results: {filepath}")
print(f"Overall verdict: {results['overall_verdict']}")
print(f"{'='*70}")
