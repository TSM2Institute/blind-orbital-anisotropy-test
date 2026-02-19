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
