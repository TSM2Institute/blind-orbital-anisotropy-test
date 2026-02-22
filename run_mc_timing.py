"""MC timing test for DESI BGS and LRG — 5 realisations each."""

import time
from boat.montecarlo import run_mc_batch

print("=" * 60)
print("MC TIMING TEST — DESI BGS (5 realisations)")
print("=" * 60)

t0 = time.time()
result_bgs = run_mc_batch("DESI_BGS", batch_size=5, total_realisations=5)
elapsed_bgs = time.time() - t0

print(f"\nDESI BGS timing test:")
print(f"  Total time for 5 realisations: {elapsed_bgs:.1f}s")

print("\n" + "=" * 60)
print("MC TIMING TEST — DESI LRG (5 realisations)")
print("=" * 60)

t0 = time.time()
result_lrg = run_mc_batch("DESI_LRG", batch_size=5, total_realisations=5)
elapsed_lrg = time.time() - t0

print(f"\nDESI LRG timing test:")
print(f"  Total time for 5 realisations: {elapsed_lrg:.1f}s")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"BGS total: {elapsed_bgs:.1f}s")
print(f"LRG total: {elapsed_lrg:.1f}s")
