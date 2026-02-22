"""
BOAT v2.0 — Survey Data Loaders and Pre-Qualification

Loads, cleans, and pre-qualifies galaxy survey data for the
blind orbital anisotropy test.

Supported surveys:
- 6dFGS DR3 (Southern, available)
- SDSS DR18 (Northern, available)
- DESI DR1 BGS (Both hemispheres, pending data download)
- DESI DR1 LRG (Both hemispheres, pending data download)

Usage:
    from boat.datasets import load_survey, get_available_surveys
    df = load_survey("6dFGS")
"""

import numpy as np
import pandas as pd
import os
import gzip

from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

from boat.manifest import (
    Z_MIN, Z_MAX, Z_MIN_LRG, Z_MAX_LRG,
    GALACTIC_PLANE_MASK_DEG, MIN_SAMPLE_SIZE,
    MIN_SKY_FRACTION, MIN_MEDIAN_Z,
    DUPLICATE_MATCH_ARCSEC, DATA_DIR, DATA_FILES,
    C_KM_S
)

LRG_SURVEYS = {"DESI_LRG"}


def load_6dfgs():
    """
    Load 6dFGS DR3 from gzipped text file.

    File format: space-separated, comment lines start with #.
    Columns (1-indexed):
        2,3,4  — RA (J2000): hours, minutes, seconds
        5,6,7  — Dec (J2000): degrees, arcminutes, arcseconds
        15     — Best combined recession velocity cz (km/s)
        16     — Uncertainty in cz (km/s)

    Returns DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, DATA_FILES["6dFGS"])

    rows = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            parts = line.split()
            if len(parts) < 16:
                continue

            ra_h = float(parts[1])
            ra_m = float(parts[2])
            ra_s = float(parts[3])
            dec_d = float(parts[4])
            dec_m = float(parts[5])
            dec_s = float(parts[6])
            cz = float(parts[14])
            cz_err = float(parts[15])

            ra_deg = 15.0 * (ra_h + ra_m / 60.0 + ra_s / 3600.0)

            dec_sign = -1.0 if dec_d < 0 else 1.0
            dec_deg = dec_sign * (abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0)

            z = cz / C_KM_S

            if cz_err > 0:
                rows.append({'ra': ra_deg, 'dec': dec_deg, 'z': z})

    df = pd.DataFrame(rows)
    return df


def load_sdss():
    """
    Load SDSS DR18 from CSV file.

    File is a CasJobs export with columns:
    ra, dec, z_obs, zErr, zWarning, class, subClass

    Returns DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, DATA_FILES["SDSS"])

    df = pd.read_csv(filepath)

    df = df.rename(columns={'z_obs': 'z'})

    df = df[['ra', 'dec', 'z']].copy()

    df = df.dropna(subset=['ra', 'dec', 'z'])

    return df


def load_desi_bgs():
    """
    Load DESI DR1 BGS survey data.

    Data file not yet available — will be downloaded from
    NOIRLab Astro Data Lab after account approval.

    Expected file: attached_assets/DESI_DR1_BGS.csv (or .fits)
    Expected columns: RA, DEC, Z (spectroscopic)
    Quality flag: ZWARN == 0

    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, "DESI_DR1_BGS.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"DESI BGS data not found at {filepath}. "
            "Download from NOIRLab Astro Data Lab first."
        )
    raise NotImplementedError("DESI BGS loader — complete after data download")


def load_desi_lrg():
    """
    Load DESI DR1 LRG survey data.

    Same source as BGS but different target class and redshift range.

    Expected file: attached_assets/DESI_DR1_LRG.csv (or .fits)

    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    filepath = os.path.join(DATA_DIR, "DESI_DR1_LRG.csv")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"DESI LRG data not found at {filepath}. "
            "Download from NOIRLab Astro Data Lab first."
        )
    raise NotImplementedError("DESI LRG loader — complete after data download")


SURVEY_LOADERS = {
    "6dFGS": load_6dfgs,
    "SDSS": load_sdss,
    "DESI_BGS": load_desi_bgs,
    "DESI_LRG": load_desi_lrg,
}


def get_available_surveys():
    """
    Return list of surveys whose data files are present.
    Does NOT check pre-qualification — just file existence.
    """
    available = []
    file_checks = {
        "6dFGS": os.path.join(DATA_DIR, DATA_FILES.get("6dFGS", "")),
        "SDSS": os.path.join(DATA_DIR, DATA_FILES.get("SDSS", "")),
        "DESI_BGS": os.path.join(DATA_DIR, "DESI_DR1_BGS.csv"),
        "DESI_LRG": os.path.join(DATA_DIR, "DESI_DR1_LRG.csv"),
    }
    for name, path in file_checks.items():
        if path and os.path.exists(path):
            available.append(name)
    return available


def load_survey(survey_name):
    """
    Load a survey by name.

    Args:
        survey_name: one of "6dFGS", "SDSS", "DESI_BGS", "DESI_LRG"

    Returns:
        DataFrame with columns: 'ra', 'dec', 'z'
    """
    if survey_name not in SURVEY_LOADERS:
        raise ValueError(f"Unknown survey: {survey_name}. "
                        f"Available: {list(SURVEY_LOADERS.keys())}")
    return SURVEY_LOADERS[survey_name]()


def apply_quality_cuts(df, survey_name):
    """
    Apply standard quality cuts to survey data.

    1. Redshift range (standard or LRG)
    2. Galactic plane mask (|b| >= 10°)
    3. Duplicate removal (match_coordinates_sky, < 1 arcsec)
    4. Minimum sample size check

    Args:
        df: DataFrame with columns 'ra', 'dec', 'z'
        survey_name: string identifying the survey

    Returns:
        DataFrame after cuts
    """
    n_start = len(df)

    if survey_name in LRG_SURVEYS:
        z_min, z_max = Z_MIN_LRG, Z_MAX_LRG
    else:
        z_min, z_max = Z_MIN, Z_MAX

    df = df[(df['z'] >= z_min) & (df['z'] <= z_max)].copy()
    n_after_z = len(df)

    coords = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg, frame='icrs')
    galactic = coords.galactic
    mask_gal = np.abs(galactic.b.deg) >= GALACTIC_PLANE_MASK_DEG
    df = df[mask_gal].copy()
    n_after_gal = len(df)

    coords_clean = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg, frame='icrs')
    idx, sep, _ = match_coordinates_sky(coords_clean, coords_clean, nthneighbor=2)
    mask_dup = sep.arcsec >= DUPLICATE_MATCH_ARCSEC
    df = df[mask_dup].copy()
    n_after_dup = len(df)

    print(f"  Quality cuts for {survey_name}:")
    print(f"    Start: {n_start:,}")
    print(f"    After z-cut ({z_min}–{z_max}): {n_after_z:,}")
    print(f"    After galactic mask (|b| >= {GALACTIC_PLANE_MASK_DEG}°): {n_after_gal:,}")
    print(f"    After duplicate removal: {n_after_dup:,}")

    return df


def check_prequalification(df, survey_name):
    """
    Check if survey meets pre-qualification criteria.

    - Sky fraction >= 15%
    - Median z >= 0.04
    - N >= 10,000

    Args:
        df: DataFrame with columns 'ra', 'dec', 'z'
        survey_name: string identifying the survey

    Returns:
        (passes: bool, sky_fraction: float, median_z: float)
    """
    ra_bins = np.floor(df['ra'].values).astype(int)
    dec_bins = np.floor(df['dec'].values + 90).astype(int)
    bin_ids = ra_bins * 180 + dec_bins
    n_occupied = len(np.unique(bin_ids))
    total_bins = 360 * 180
    sky_fraction = n_occupied / total_bins

    median_z = np.median(df['z'].values)

    passes = (sky_fraction >= MIN_SKY_FRACTION) and (median_z >= MIN_MEDIAN_Z) and (len(df) >= MIN_SAMPLE_SIZE)

    print(f"  Pre-qualification for {survey_name}:")
    print(f"    Sky fraction: {sky_fraction:.1%} (threshold: {MIN_SKY_FRACTION:.0%}) — {'PASS' if sky_fraction >= MIN_SKY_FRACTION else 'FAIL'}")
    print(f"    Median z: {median_z:.4f} (threshold: {MIN_MEDIAN_Z}) — {'PASS' if median_z >= MIN_MEDIAN_Z else 'FAIL'}")
    print(f"    Sample size: {len(df):,} (threshold: {MIN_SAMPLE_SIZE:,}) — {'PASS' if len(df) >= MIN_SAMPLE_SIZE else 'FAIL'}")
    print(f"    Overall: {'QUALIFIED' if passes else 'DISQUALIFIED'}")

    return passes, sky_fraction, median_z
