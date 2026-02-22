# BOAT v2.0 — DESI Loaders + MC Timing Test

**Context:** DESI DR1 BGS (5.9M galaxies) and LRG (4.6M galaxies) downloaded as clean CSV with columns `ra, dec, z`. All quality filters applied server-side. Files in `attached_assets/`.

---

## TASK 1: Update DESI loaders in datasets.py

Replace the `load_desi_bgs()` and `load_desi_lrg()` stubs with working implementations. These are simple — the files are already clean CSV with our standard column names.

```python
def load_desi_bgs():
    """Load DESI DR1 BGS survey data (bright programme)."""
    filepath = os.path.join(DATA_DIR, "DESI_DR1_BGS.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"DESI BGS data not found at {filepath}")
    df = pd.read_csv(filepath)
    # Columns are already 'ra', 'dec', 'z' — no renaming needed
    return df


def load_desi_lrg():
    """Load DESI DR1 LRG survey data (dark programme)."""
    filepath = os.path.join(DATA_DIR, "DESI_DR1_LRG.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"DESI LRG data not found at {filepath}")
    df = pd.read_csv(filepath)
    return df
```

Also update the `DATA_FILES` dict in `boat/manifest.py` to include the DESI filenames:

```python
DATA_FILES = {
    "6dFGS": "6dFGS_DR3.txt.gz",
    "SDSS": "SDSS_DR18.csv",
    "DESI_BGS": "DESI_DR1_BGS.csv",
    "DESI_LRG": "DESI_DR1_LRG.csv",
}
```

And update the `file_checks` dict in `get_available_surveys()` in datasets.py to use the manifest DATA_FILES for DESI too:

```python
def get_available_surveys():
    available = []
    file_checks = {
        "6dFGS": os.path.join(DATA_DIR, DATA_FILES.get("6dFGS", "")),
        "SDSS": os.path.join(DATA_DIR, DATA_FILES.get("SDSS", "")),
        "DESI_BGS": os.path.join(DATA_DIR, DATA_FILES.get("DESI_BGS", "")),
        "DESI_LRG": os.path.join(DATA_DIR, DATA_FILES.get("DESI_LRG", "")),
    }
    for name, path in file_checks.items():
        if path and os.path.exists(path):
            available.append(name)
    return available
```

Also make sure `SURVEY_LOADERS` includes the DESI entries (not commented out):

```python
SURVEY_LOADERS = {
    "6dFGS": load_6dfgs,
    "SDSS": load_sdss,
    "DESI_BGS": load_desi_bgs,
    "DESI_LRG": load_desi_lrg,
}
```

---

## TASK 2: Verify DESI data loads and pre-qualifies

```python
from boat.datasets import load_survey, apply_quality_cuts, check_prequalification, get_available_surveys

print(f"Available surveys: {get_available_surveys()}")

for name in ["DESI_BGS", "DESI_LRG"]:
    print(f"\n{'='*50}")
    print(f"Loading {name}...")
    df = load_survey(name)
    print(f"  Raw: {len(df):,} galaxies")
    
    df = apply_quality_cuts(df, name)
    print(f"  After cuts: {len(df):,}")
    
    passes, sky_frac, med_z = check_prequalification(df, name)
    print(f"  Pre-qualified: {passes}")
```

Report the galaxy counts after quality cuts, sky fractions, and median z for both DESI surveys.

---

## TASK 3: MC Timing Test (CRITICAL)

Before committing to 10,000 realisations on 5.9M galaxies, we need actual timing data.

Run a **5-realisation smoke test** on DESI BGS to measure:
- Footprint weight computation time
- Per-realisation time
- Memory usage (approximate)

```python
import time
from boat.montecarlo import run_mc_batch

# Just 5 realisations to get timing
t0 = time.time()
result = run_mc_batch("DESI_BGS", batch_size=5, total_realisations=5)
elapsed = time.time() - t0

print(f"\nDESI BGS timing test:")
print(f"  Total time for 5 realisations: {elapsed:.1f}s")
print(f"  Estimated per-realisation: {(elapsed - 238) / 5:.1f}s")  # subtract approx weight time
print(f"  Estimated time for 1,000 realisations: {((elapsed - 238) / 5 * 1000 + 238) / 60:.0f} minutes")
print(f"  Estimated time for 10,000 realisations: {((elapsed - 238) / 5 * 10000 + 238) / 3600:.1f} hours")
```

Also run on DESI LRG:

```python
t0 = time.time()
result_lrg = run_mc_batch("DESI_LRG", batch_size=5, total_realisations=5)
elapsed_lrg = time.time() - t0

print(f"\nDESI LRG timing test:")
print(f"  Total time for 5 realisations: {elapsed_lrg:.1f}s")
```

---

## REPORT BACK

1. DESI loaders updated: yes/no
2. manifest.py DATA_FILES updated: yes/no
3. get_available_surveys() now returns all 4: yes/no
4. DESI BGS after quality cuts: galaxy count, sky fraction, median z
5. DESI LRG after quality cuts: galaxy count, sky fraction, median z
6. Both pre-qualify: yes/no
7. **BGS timing test:**
   - Weight computation time
   - Per-realisation time
   - Estimated time for 10,000 realisations
8. **LRG timing test:**
   - Weight computation time
   - Per-realisation time
   - Estimated time for 10,000 realisations

**Based on timing results, we'll decide whether to run full 10k, reduce realisations, or subsample galaxies for MC calibration.**

Wait for instruction before starting any long MC runs.

---

*Build prompt prepared by Claude Opus (Anthropic), 22 February 2026*
