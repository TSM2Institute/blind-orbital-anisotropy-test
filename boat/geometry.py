"""
BOAT Geometry Module
Converts RA/Dec → unit vector → velocity projection → SR Doppler z_model.
No observed redshift is used at any point in this module.
"""

import numpy as np
from boat.manifest import V_OBS_VEC, C_KM_S


def radec_to_unit_vector(ra_deg, dec_deg):
    """
    Convert RA/Dec (degrees) to Cartesian unit vector(s) in ICRS.

    Parameters
    ----------
    ra_deg : float or array
        Right Ascension in degrees
    dec_deg : float or array
        Declination in degrees

    Returns
    -------
    numpy.ndarray
        Unit vector(s), shape (3,) for scalar input or (N, 3) for array input.
    """
    ra = np.radians(np.atleast_1d(ra_deg))
    dec = np.radians(np.atleast_1d(dec_deg))

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    vectors = np.column_stack([x, y, z])

    if vectors.shape[0] == 1:
        return vectors[0]
    return vectors


def compute_v_los(n_hat):
    """
    Compute line-of-sight velocity for galaxy direction(s).

    v_LOS = -v_obs · n_hat
    Negative sign: motion toward source = positive v_LOS = blueshift.

    Parameters
    ----------
    n_hat : numpy.ndarray
        Unit vector(s), shape (3,) or (N, 3)

    Returns
    -------
    float or numpy.ndarray
        Line-of-sight velocity in km/s. Positive = recession (redshift),
        negative = approach (blueshift).
    """
    if n_hat.ndim == 1:
        return -np.dot(V_OBS_VEC, n_hat)
    return -np.dot(n_hat, V_OBS_VEC)


def v_los_to_z_model(v_los):
    """
    Convert line-of-sight velocity to SR Doppler redshift.

    z_model = sqrt((1 + beta) / (1 - beta)) - 1
    where beta = v_los / c

    Parameters
    ----------
    v_los : float or numpy.ndarray
        Line-of-sight velocity in km/s

    Returns
    -------
    float or numpy.ndarray
        Predicted redshift (positive = redshift, negative = blueshift)
    """
    beta = v_los / C_KM_S
    beta = np.clip(beta, -0.9999, 0.9999)
    return np.sqrt((1 + beta) / (1 - beta)) - 1


def predict_z_model(ra_deg, dec_deg):
    """
    Full pipeline: RA/Dec → z_model (no observed redshift used).

    Parameters
    ----------
    ra_deg : float or array
        Right Ascension in degrees
    dec_deg : float or array
        Declination in degrees

    Returns
    -------
    float or numpy.ndarray
        Predicted z_model from orbital geometry alone
    """
    n_hat = radec_to_unit_vector(ra_deg, dec_deg)
    v_los = compute_v_los(n_hat)
    return v_los_to_z_model(v_los)


def angular_separation(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """
    Compute angular separation between two sky positions (degrees).
    Uses the Vincenty formula for numerical stability.

    Parameters
    ----------
    ra1_deg, dec1_deg : float
        First position (degrees)
    ra2_deg, dec2_deg : float
        Second position (degrees)

    Returns
    -------
    float
        Angular separation in degrees
    """
    ra1, dec1 = np.radians(ra1_deg), np.radians(dec1_deg)
    ra2, dec2 = np.radians(ra2_deg), np.radians(dec2_deg)

    dra = ra2 - ra1

    num = np.sqrt(
        (np.cos(dec2) * np.sin(dra))**2 +
        (np.cos(dec1) * np.sin(dec2) - np.sin(dec1) * np.cos(dec2) * np.cos(dra))**2
    )
    den = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(dra)

    return np.degrees(np.arctan2(num, den))
