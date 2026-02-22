"""Run SDSS full analysis with z-bin and top-10 axes diagnostics."""

from boat.montecarlo import load_sdss, apply_quality_cuts
from boat.analysis import run_full_analysis

df = load_sdss()
df = apply_quality_cuts(df, "SDSS")

results = run_full_analysis(
    df['ra'].values, df['dec'].values, df['z'].values, "SDSS"
)

print("\nz-bin invariance details:")
for b in results['zbin_invariance']['bins']:
    if 'r' in b and b['r'] is not None:
        print(f"  {b['z_range']}: N={b['n_galaxies']:,}, r={b['r']:.6f}, |r|={b['abs_r']:.6f}, p={b['p_value']:.2e}")
    else:
        print(f"  {b['z_range']}: N={b['n_galaxies']:,}, {b['status']}")

print("\nTop 10 axes by |r|:")
sorted_axes = sorted(results['all_axes'], key=lambda x: x['abs_r'], reverse=True)
for i, a in enumerate(sorted_axes[:10]):
    print(f"  #{i+1}: {a['label']:20s} |r|={a['abs_r']:.6f} p={a['p_value']:.2e}")

print(f"\nSDSS verdict: {results['verdict']}")
print(f"Ratio: {results['primary_test']['ratio']:.4f}")
print(f"Rank: #{results['primary_test']['rank']}")
print(f"Orbital discriminator: {results['orbital_discriminator']['verdict']}")
