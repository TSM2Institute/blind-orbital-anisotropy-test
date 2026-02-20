"""
BOAT Dataset Module
Loads, filters, and prepares galaxy survey datasets.
Applies all quality cuts from the sealed manifest.
"""

import numpy as np
import pandas as pd
import os
from astropy.coordinates import SkyCoord
import astropy.units as u

from boat.manifest import (
    Z_MIN, Z_MAX, GALACTIC_PLANE_MASK_DEG,
    MIN_SAMPLE_SIZE, DUPLICATE_MATCH_ARCSEC,
    DATASET_SELECTION_SEED, N_DATASETS
)


SURVEY_REGISTRY = {
    "synthetic_uniform": {
        "description": "Synthetic uniform sky coverage for testing",
        "type": "synthetic",
        "n_objects": 50000,
    },
    "synthetic_dipole": {
        "description": "Synthetic data WITH injected dipole signal for validation",
        "type": "synthetic_dipole",
        "n_objects": 50000,
        "injected_amplitude": 0.00124,
    },
    "6dFGS_DR3": {
        "description": "6dF Galaxy Survey Data Release 3 (final)",
        "type": "real",
        "filepath": "attached_assets/6dFGSzDR3.txt_1771500188269.gz",
        "n_expected": 124647,
    },
    "2dFGRS": {
        "description": "2dF Galaxy Redshift Survey (final release)",
        "type": "real",
        "filepath": "attached_assets/best.observations.idz_1771582454022.gz",
        "n_expected": 245591,
    },
    "GAMA_DR4": {
        "description": "Galaxy and Mass Assembly DR4 SpecObj",
        "type": "real",
        "filepath": "attached_assets/SpecObjv27_1771583256206.fits",
        "n_expected": 344905,
    },
    "2MRS": {
        "description": "2MASS Redshift Survey (Huchra+ 2012)",
        "type": "real",
        "filepath": "attached_assets/asu_1771586798714.tsv",
        "n_expected": 44599,
    },
}


def generate_synthetic_uniform(n_objects, seed):
    """
    Generate synthetic galaxy catalogue with uniform sky coverage.
    No dipole signal injected — useful for null testing.

    Parameters
    ----------
    n_objects : int
        Number of galaxies to generate
    seed : int
        Random seed for reproducibility

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, source
    """
    rng = np.random.RandomState(seed)

    ra = rng.uniform(0, 360, n_objects)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n_objects)))

    z_obs = rng.uniform(Z_MIN, Z_MAX, n_objects)

    return pd.DataFrame({
        'ra': ra,
        'dec': dec,
        'z_obs': z_obs,
        'source': 'synthetic_uniform'
    })


def generate_synthetic_dipole(n_objects, seed, amplitude=0.00124):
    """
    Generate synthetic galaxy catalogue with an INJECTED dipole signal
    matching the expected BOAT orbital signature. Used to validate
    that the analysis pipeline can detect a known signal.

    Parameters
    ----------
    n_objects : int
        Number of galaxies
    seed : int
        Random seed
    amplitude : float
        Dipole amplitude to inject (in z units)

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, z_injected, source
    """
    from boat.geometry import predict_z_model

    rng = np.random.RandomState(seed)

    ra = rng.uniform(0, 360, n_objects)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n_objects)))

    z_base = rng.uniform(Z_MIN, Z_MAX, n_objects)

    z_dipole = predict_z_model(ra, dec)

    scale_factor = amplitude / 0.00124
    z_injected = z_dipole * scale_factor

    z_noise = rng.normal(0, 0.001, n_objects)

    z_obs = z_base + z_injected + z_noise

    return pd.DataFrame({
        'ra': ra,
        'dec': dec,
        'z_obs': z_obs,
        'z_injected': z_injected,
        'source': 'synthetic_dipole'
    })


def angular_separation_simple(ra1, dec1, ra2, dec2):
    """Quick angular separation in degrees (scalar only)."""
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)
    cos_sep = (np.sin(dec1) * np.sin(dec2) +
               np.cos(dec1) * np.cos(dec2) * np.cos(ra2 - ra1))
    return np.degrees(np.arccos(np.clip(cos_sep, -1, 1)))


def apply_quality_cuts(df, source_name="unknown"):
    """
    Apply all manifest-specified quality cuts to a dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: ra, dec, z_obs
    source_name : str
        Dataset name for logging

    Returns
    -------
    pd.DataFrame
        Filtered dataset
    dict
        Filter statistics
    """
    stats = {'source': source_name, 'n_raw': len(df)}

    mask = (
        np.isfinite(df['ra']) & np.isfinite(df['dec']) &
        (df['ra'] >= 0) & (df['ra'] < 360) &
        (df['dec'] >= -90) & (df['dec'] <= 90)
    )
    df = df[mask].copy()
    stats['n_after_coord_cut'] = len(df)

    mask = (df['z_obs'] >= Z_MIN) & (df['z_obs'] <= Z_MAX)
    df = df[mask].copy()
    stats['n_after_z_cut'] = len(df)

    coords = SkyCoord(ra=df['ra'].values * u.deg,
                      dec=df['dec'].values * u.deg, frame='icrs')
    gal_b = coords.galactic.b.deg
    mask = np.abs(gal_b) >= GALACTIC_PLANE_MASK_DEG
    df = df[mask].copy()
    df['gal_b'] = gal_b[mask]
    stats['n_after_gal_cut'] = len(df)

    if len(df) > 1:
        keep = np.ones(len(df), dtype=bool)
        ra_sorted_idx = np.argsort(df['ra'].values)
        ra_vals = df['ra'].values[ra_sorted_idx]
        dec_vals = df['dec'].values[ra_sorted_idx]

        for i in range(1, len(ra_sorted_idx)):
            if not keep[ra_sorted_idx[i]]:
                continue
            if abs(ra_vals[i] - ra_vals[i-1]) > DUPLICATE_MATCH_ARCSEC / 3600 * 2:
                continue
            sep = angular_separation_simple(
                ra_vals[i], dec_vals[i], ra_vals[i-1], dec_vals[i-1]
            )
            if sep < DUPLICATE_MATCH_ARCSEC / 3600:
                keep[ra_sorted_idx[i]] = False

        df = df[keep].copy()
    stats['n_after_dedup'] = len(df)

    stats['passes_min_size'] = len(df) >= MIN_SAMPLE_SIZE
    stats['n_final'] = len(df)

    return df, stats


def load_6dfgs(filepath='attached_assets/6dFGSzDR3.txt_1771500188269.gz'):
    """
    Load 6dF Galaxy Survey DR3 redshift catalogue.

    File format: gzipped space-separated ASCII, 75 header lines starting with #.
    Key columns (0-indexed):
        0: target ID
        1,2,3: RA hours, minutes, seconds (J2000)
        4,5,6: Dec degrees, arcminutes, arcseconds (J2000)
        14: best cz (km/s)
        15: cz uncertainty (km/s)
        17: redshift quality flag (3 or 4 = reliable for 6dFGS sources)
        16: source code (126 = 6dFGS)
        18: galactic latitude
        19: galactic longitude

    Parameters
    ----------
    filepath : str
        Path to 6dFGSzDR3.txt.gz

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, source
    """
    import gzip

    rows = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 20:
                continue

            try:
                ra_h = float(parts[1])
                ra_m = float(parts[2])
                ra_s = float(parts[3])
                ra_deg = (ra_h + ra_m / 60 + ra_s / 3600) * 15.0

                dec_d = float(parts[4])
                dec_m = float(parts[5])
                dec_s = float(parts[6])
                sign = -1 if dec_d < 0 or parts[4].startswith('-') else 1
                dec_deg = sign * (abs(dec_d) + dec_m / 60 + dec_s / 3600)

                cz = float(parts[14])
                cz_err = float(parts[15])
                source_code = int(parts[16])
                quality = int(parts[17])

                z_obs = cz / 299792.458

                rows.append({
                    'ra': ra_deg,
                    'dec': dec_deg,
                    'z_obs': z_obs,
                    'cz_err': cz_err,
                    'quality': quality,
                    'source_code': source_code,
                    'source': '6dFGS_DR3'
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)

    mask_6df = (df['source_code'] == 126) & (df['quality'].isin([3, 4]))
    mask_other = df['source_code'] != 126
    df = df[mask_6df | mask_other].copy()

    df = df[df['z_obs'] > 0].copy()

    print(f"  6dFGS DR3 loaded: {len(df)} galaxies with reliable redshifts")
    print(f"  z range: [{df['z_obs'].min():.4f}, {df['z_obs'].max():.4f}]")
    print(f"  Median z: {df['z_obs'].median():.4f}")

    return df[['ra', 'dec', 'z_obs', 'source']].copy()


def load_2dfgrs(filepath='attached_assets/best.observations.idz_1771582454022.gz'):
    """
    Load 2dF Galaxy Redshift Survey final catalogue.

    File format: gzipped ASCII, 35 whitespace-delimited fields.
    Key columns (0-indexed):
        4,5,6: RA J2000 (h, m, s)
        7,8,9: Dec J2000 (d, arcmin, arcsec)
        28: Best redshift
        29: Quality of best redshift (>=3 reliable)

    Sentinel values: -9.999x = missing data

    Parameters
    ----------
    filepath : str
        Path to best.observations.idz.gz

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, source
    """
    import gzip

    rows = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 30:
                continue

            try:
                ra_h = float(parts[4])
                ra_m = float(parts[5])
                ra_s = float(parts[6])
                ra_deg = (ra_h + ra_m / 60 + ra_s / 3600) * 15.0

                dec_d = float(parts[7])
                dec_m = float(parts[8])
                dec_s = float(parts[9])
                sign = -1 if dec_d < 0 or parts[7].startswith('-') else 1
                dec_deg = sign * (abs(dec_d) + dec_m / 60 + dec_s / 3600)

                z_best = float(parts[28])
                quality = float(parts[29])

                if z_best < -9 or quality < 0:
                    continue

                if quality < 3:
                    continue

                if z_best <= 0:
                    continue

                rows.append({
                    'ra': ra_deg,
                    'dec': dec_deg,
                    'z_obs': z_best,
                    'quality': int(quality),
                    'source': '2dFGRS'
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)

    print(f"  2dFGRS loaded: {len(df)} galaxies with reliable redshifts (quality >= 3)")
    print(f"  z range: [{df['z_obs'].min():.4f}, {df['z_obs'].max():.4f}]")
    print(f"  Median z: {df['z_obs'].median():.4f}")

    return df[['ra', 'dec', 'z_obs', 'source']].copy()


def load_2mrs(filepath='attached_assets/asu_1771586798714.tsv'):
    """
    Load 2MASS Redshift Survey (2MRS) from VizieR TSV export.

    Parses table3 only (main catalogue, ~44,599 galaxies).
    Tab-separated, column headers at line 66, data from line 69.

    Key columns: RAJ2000 (deg), DEJ2000 (deg), cz (km/s heliocentric)
    z_obs = cz / c

    Parameters
    ----------
    filepath : str
        Path to the VizieR TSV file

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, source
    """
    rows = []
    in_table3 = False
    headers = None
    past_separator = False

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()

            if 'J/ApJS/199/26/table3' in stripped:
                in_table3 = True
                past_separator = False
                continue

            if in_table3 and stripped.startswith('#') and 'J/ApJS/199/26/' in stripped and 'table3' not in stripped:
                break

            if not in_table3:
                continue

            if stripped.startswith('#'):
                continue

            if headers is None and not stripped.startswith('-') and stripped:
                headers = stripped.split('\t')
                headers = [h.strip() for h in headers]
                continue

            if not past_separator:
                if stripped.startswith('-'):
                    past_separator = True
                continue

            if not stripped:
                continue

            parts = stripped.split('\t')

            try:
                ra_idx = headers.index('RAJ2000') if 'RAJ2000' in headers else None
                dec_idx = headers.index('DEJ2000') if 'DEJ2000' in headers else None
                cz_idx = headers.index('cz') if 'cz' in headers else None

                if ra_idx is None or dec_idx is None or cz_idx is None:
                    continue

                ra = float(parts[ra_idx].strip())
                dec = float(parts[dec_idx].strip())
                cz = float(parts[cz_idx].strip())

                z_obs = cz / 299792.458

                rows.append({
                    'ra': ra,
                    'dec': dec,
                    'z_obs': z_obs,
                    'source': '2MRS'
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)

    df = df[np.isfinite(df['z_obs']) & np.isfinite(df['ra']) & np.isfinite(df['dec'])].copy()

    print(f"  2MRS loaded: {len(df)} galaxies")
    print(f"  z range: [{df['z_obs'].min():.4f}, {df['z_obs'].max():.4f}]")
    print(f"  Median z: {df['z_obs'].median():.4f}")

    return df[['ra', 'dec', 'z_obs', 'source']].copy()


def load_gama(filepath='attached_assets/SpecObjv27_1771583256206.fits'):
    """
    Load GAMA DR4 SpecObj catalogue.

    FITS binary table. Columns: RA, DEC (J2000 degrees), Z (redshift),
    NQ (quality flag, >=3 reliable).

    Parameters
    ----------
    filepath : str
        Path to SpecObj.fits

    Returns
    -------
    pd.DataFrame
        Columns: ra, dec, z_obs, source
    """
    from astropy.io import fits
    from astropy.table import Table

    t = Table.read(filepath, hdu='SpecObj')
    df = t.to_pandas()

    print(f"  GAMA DR4 raw: {len(df)} rows")

    df = df[df['NQ'] >= 3].copy()

    df = df[df['Z'] > 0].copy()
    df = df[np.isfinite(df['Z'])].copy()

    df = df[np.isfinite(df['RA']) & np.isfinite(df['DEC'])].copy()

    df = df.rename(columns={'RA': 'ra', 'DEC': 'dec', 'Z': 'z_obs'})
    df['source'] = 'GAMA_DR4'

    print(f"  GAMA DR4 after quality cuts: {len(df)} galaxies (NQ >= 3)")
    print(f"  z range: [{df['z_obs'].min():.4f}, {df['z_obs'].max():.4f}]")
    print(f"  Median z: {df['z_obs'].median():.4f}")

    return df[['ra', 'dec', 'z_obs', 'source']].copy()


def load_dataset(name, seed=None):
    """
    Load a dataset by name from the registry.

    Parameters
    ----------
    name : str
        Dataset name (must be in SURVEY_REGISTRY)
    seed : int, optional
        Random seed (for synthetic datasets)

    Returns
    -------
    pd.DataFrame
        Raw dataset before quality cuts
    """
    if name not in SURVEY_REGISTRY:
        raise ValueError(f"Unknown dataset: {name}. Available: {list(SURVEY_REGISTRY.keys())}")

    info = SURVEY_REGISTRY[name]

    if info['type'] == 'synthetic':
        return generate_synthetic_uniform(info['n_objects'], seed or 42)
    elif info['type'] == 'synthetic_dipole':
        return generate_synthetic_dipole(
            info['n_objects'], seed or 42,
            amplitude=info.get('injected_amplitude', 0.00124)
        )
    elif info['type'] == 'real':
        filepath = info.get('filepath', '')
        if '6dFGS' in name or '6dfgs' in name.lower():
            return load_6dfgs(filepath)
        elif '2dFGRS' in name or '2dfgrs' in name.lower():
            return load_2dfgrs(filepath)
        elif 'GAMA' in name or 'gama' in name.lower():
            return load_gama(filepath)
        elif '2MRS' in name or '2mrs' in name.lower():
            return load_2mrs(filepath)
        else:
            raise NotImplementedError(f"No loader for real dataset '{name}'")
    else:
        raise NotImplementedError(f"Dataset type '{info['type']}' not yet implemented")


def select_datasets(seed=DATASET_SELECTION_SEED):
    """
    Randomly select K datasets from the registry.
    For development: returns all available datasets.
    For production: will randomly select from real surveys.

    Parameters
    ----------
    seed : int
        Random seed for reproducible selection

    Returns
    -------
    list of str
        Selected dataset names
    """
    available = list(SURVEY_REGISTRY.keys())

    if len(available) <= N_DATASETS:
        return available

    rng = np.random.RandomState(seed)
    return list(rng.choice(available, size=N_DATASETS, replace=False))
