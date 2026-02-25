# BOAT v2.0 — Manifest Sealing Prompt

**Context:** R_NULL_99 = 1.541144 is locked in manifest.py and verified. 40,000 MC realisations complete across 4 surveys. All pipeline modules built and tested. This prompt seals the manifest document.

**THIS IS A CRITICAL STEP. Execute exactly as written. No creative additions.**

---

## TASK: Update and seal BOAT_v2_MANIFEST_FINAL_DRAFT.md

The manifest document needs to be uploaded to the Replit project if it isn't already there. Graham will provide it. Once it's in the project, make exactly these four edits:

### Edit 1: Date field

Find:
```
**Date:** [TO BE SEALED AFTER MONTE CARLO CALIBRATION]
```

Replace with:
```
**Date:** 23 February 2026 (sealed)
```

### Edit 2: DESI access column in Section 5.2

Find:
```
| 3 | DESI DR1 BGS | ~18% | ~0.20 | Both (northern-heavy) | NERSC |
| 4 | DESI DR1 LRG | ~18% | ~0.70 | Both (northern-heavy) | NERSC |
```

Replace with:
```
| 3 | DESI DR1 BGS | ~18% | ~0.20 | Both (northern-heavy) | NOIRLab Astro Data Lab |
| 4 | DESI DR1 LRG | ~18% | ~0.70 | Both (northern-heavy) | NOIRLab Astro Data Lab |
```

### Edit 3: Section 7.8 — Fill R_NULL_99 and add subsampling note

Find:
```
**Locked value:** R_NULL_99 = [TO BE FILLED BY MONTE CARLO — locked before
any real v2 analysis]

**Documentation:** The full distribution of null ratios, per-survey
breakdowns, and the 95th/99th/99.5th percentile values are recorded
in the Monte Carlo results file (SHA-256 hashed and committed).
```

Replace with:
```
**Locked value:** R_NULL_99 = 1.541144

Determined from 40,000 combined null realisations (10,000 per survey × 4 surveys).

Per-survey 99th percentiles:

| Survey | N galaxies (MC) | Realisations | p99 |
|--------|----------------|--------------|------|
| 6dFGS DR3 | 119,906 | 10,000 | 1.5776 |
| SDSS DR18 | 1,575,204 | 10,000 | 1.4984 |
| DESI DR1 BGS | 1,000,000 (subsample) | 10,000 | 1.5369 |
| DESI DR1 LRG | 1,000,000 (subsample) | 10,000 | 1.5264 |
| **Combined** | — | **40,000** | **1.5411** |

For DESI BGS and DESI LRG, Monte Carlo null realisations used a
1,000,000-galaxy random subsample (seed 12345) of the quality-cut
catalogue, preserving the survey footprint geometry. Validation
confirmed the per-survey p99 values are consistent with full-catalogue
calibration on 6dFGS and SDSS (range: 1.498–1.578 across all four
surveys). Production analysis uses the complete catalogues.

**Documentation:** The full distribution of null ratios, per-survey
breakdowns, and the 95th/99th/99.5th percentile values are recorded
in the Monte Carlo results file (SHA-256 hashed and committed).
```

### Edit 4: Section 13 — Compute and insert SHA-256 hash

After making Edits 1–3, save the file as `BOAT_v2_MANIFEST.md` (drop the `_FINAL_DRAFT` suffix — it's now sealed).

Then compute the SHA-256:

```python
import hashlib

with open("BOAT_v2_MANIFEST.md", "rb") as f:
    manifest_hash = hashlib.sha256(f.read()).hexdigest()

print(f"Manifest SHA-256: {manifest_hash}")
```

Now edit the manifest one final time — find:

```
**Manifest hash (pre-seal):** [COMPUTED AT SEAL TIME, AFTER R_NULL_99 IS FILLED]
```

Replace with:

```
**Manifest hash (pre-seal):** [SHA-256 hash computed from this document with this field blank]
```

**IMPORTANT:** The hash is computed on the document BEFORE this field is filled. This is standard practice — the hash covers everything except itself. Record the hash value in the results JSON and in the GitHub commit message.

Actually, simpler approach: compute the hash with the placeholder text still in place, then record it separately. The hash authenticates the full document content.

Use this procedure instead:

1. Make Edits 1, 2, 3
2. Leave Edit 4 placeholder as-is: `[COMPUTED AT SEAL TIME, AFTER R_NULL_99 IS FILLED]`
3. Save as `BOAT_v2_MANIFEST.md`
4. Compute SHA-256 of the complete file
5. Print the hash — Graham will record it in the GitHub commit

```python
import hashlib

filepath = "BOAT_v2_MANIFEST.md"

with open(filepath, "rb") as f:
    content = f.read()
    manifest_hash = hashlib.sha256(content).hexdigest()

print(f"\nManifest file: {filepath}")
print(f"Size: {len(content):,} bytes")
print(f"SHA-256: {manifest_hash}")
print(f"\nThis hash authenticates the sealed manifest.")
print(f"Record in GitHub commit message and production results JSON.")
```

---

## REPORT BACK

1. Edit 1 (date): done
2. Edit 2 (DESI access): done
3. Edit 3 (R_NULL_99 + subsampling note): done
4. File renamed to BOAT_v2_MANIFEST.md: done
5. **SHA-256 hash of sealed manifest:** [the hash]

**Do not run production yet.** Wait for Graham to confirm the sealed manifest before we execute.

---

*Sealing prompt prepared by Claude Opus (Anthropic), 22 February 2026*
