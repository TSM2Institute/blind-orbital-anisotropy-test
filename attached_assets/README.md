# BOAT Survey Data Files

These data files are too large for GitHub and must be downloaded separately.
Place them in this directory (`attached_assets/`) with the filenames shown below.

## Required Catalogues

| File | Size | Source | Download |
|------|------|--------|----------|
| `6dFGSzDR3.txt_1771500188269.gz` | 5.3 MB | 6dF Galaxy Survey DR3 | [VizieR: VII/259](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=VII/259/6dfgs3) |
| `best.observations.idz_1771582454022.gz` | 14 MB | 2dF Galaxy Redshift Survey | [2dFGRS Archive](http://www.2dfgrs.net/) |
| `SpecObjv27_1771583256206.fits` | 120 MB | GAMA DR4 SpecObj | [GAMA DR4](http://www.gama-survey.org/dr4/) |
| `asu_1771586798714.tsv` | 13 MB | 2MASS Redshift Survey | [VizieR: J/ApJS/199/26](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/ApJS/199/26/table3) |
| `MyTable_TSM2_0_1771642595705.csv` | 88 MB | SDSS DR18 Spectroscopic Galaxies | [SDSS CasJobs](https://skyserver.sdss.org/CasJobs/) |

## SDSS CasJobs Query

The SDSS file was generated with this SQL query on CasJobs (context: DR18):

```sql
SELECT s.ra, s.dec, s.z AS z_obs, s.zErr, s.zWarning, s.class, s.subClass
FROM SpecObj AS s
WHERE s.class = 'GALAXY'
  AND s.zWarning = 0
  AND s.z BETWEEN 0.01 AND 0.5
  AND s.zErr < 0.001
```

## Notes

- All coordinates are J2000 equatorial (RA/Dec in degrees)
- The BOAT loaders in `boat/datasets.py` expect these exact filenames
- Total data volume: ~241 MB
