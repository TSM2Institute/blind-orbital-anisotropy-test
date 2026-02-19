"""
BOAT Analysis Module
Residual computation, dipole+quadrupole regression, decoy axis comparison.
Implements the complete statistical test specification from the sealed manifest.
"""

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from boat.manifest import (
    V_HAT, N_DECOY_AXES, DECOY_AXES_SEED,
    CORRELATION_THRESHOLD, P_VALUE_THRESHOLD
)
from boat.geometry import radec_to_unit_vector, predict_z_model


def compute_residuals(df):
    """
    Compute redshift residuals by removing the isotropic component.

    delta_z_i = z_obs_i - mean(z_obs)

    Subtracting the sample mean removes the dominant Hubble-flow component,
    isolating any directional (dipole) signal in the data.

    z_model is computed for reference/comparison but is NOT subtracted
    from the data — it represents the predicted signal we're testing for.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: ra, dec, z_obs

    Returns
    -------
    pd.DataFrame
        Input dataframe with added columns: z_model, delta_z
    """
    df = df.copy()

    df['z_model'] = predict_z_model(df['ra'].values, df['dec'].values)

    z_mean = df['z_obs'].mean()
    df['delta_z'] = df['z_obs'] - z_mean

    return df


def compute_cos_theta(ra, dec, axis_vector):
    """
    Compute cos(theta) between galaxy directions and a given axis.

    Parameters
    ----------
    ra, dec : array-like
        Galaxy positions in degrees
    axis_vector : numpy.ndarray
        Unit vector defining the axis, shape (3,)

    Returns
    -------
    numpy.ndarray
        cos(theta) values for each galaxy
    """
    n_hat = radec_to_unit_vector(ra, dec)
    return np.dot(n_hat, axis_vector)


def fit_dipole(cos_theta, delta_z):
    """
    Fit dipole model: delta_z = A * cos(theta) + intercept

    Parameters
    ----------
    cos_theta : numpy.ndarray
        cos(angular separation) from axis
    delta_z : numpy.ndarray
        Redshift residuals

    Returns
    -------
    dict
        'A': dipole amplitude
        'r': Pearson correlation coefficient
        'p_value': p-value for correlation
        'rms': RMS of residuals after fit
    """
    r, p_value = sp_stats.pearsonr(cos_theta, delta_z)

    slope, intercept, _, _, std_err = sp_stats.linregress(cos_theta, delta_z)

    fitted = slope * cos_theta + intercept
    rms = np.sqrt(np.mean((delta_z - fitted) ** 2))

    return {
        'A_dipole': slope,
        'intercept': intercept,
        'r': r,
        'p_value': p_value,
        'rms': rms,
        'std_err': std_err
    }


def fit_dipole_quadrupole(cos_theta, delta_z):
    """
    Fit dipole + quadrupole model:
    delta_z = A * cos(theta) + B * (3*cos^2(theta) - 1)/2 + intercept

    Parameters
    ----------
    cos_theta : numpy.ndarray
        cos(angular separation) from axis
    delta_z : numpy.ndarray
        Redshift residuals

    Returns
    -------
    dict
        Fit results including amplitudes, significances, and diagnostics
    """
    P2 = (3 * cos_theta**2 - 1) / 2
    X = np.column_stack([cos_theta, P2, np.ones_like(cos_theta)])

    XtX = X.T @ X
    Xty = X.T @ delta_z

    try:
        beta = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        return {
            'A_dipole': np.nan, 'B_quadrupole': np.nan,
            'intercept': np.nan, 'r_squared': np.nan,
            'p_dipole': 1.0, 'p_quadrupole': 1.0,
            'rms': np.nan, 'significant_terms': []
        }

    A_dipole = beta[0]
    B_quadrupole = beta[1]
    intercept = beta[2]

    y_pred = X @ beta
    residuals = delta_z - y_pred
    rms = np.sqrt(np.mean(residuals ** 2))

    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((delta_z - np.mean(delta_z)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    n = len(delta_z)
    p = X.shape[1]
    if n > p:
        mse = ss_res / (n - p)
        try:
            cov = mse * np.linalg.inv(XtX)
            se = np.sqrt(np.diag(cov))
            t_stats = beta / se
            p_values = 2 * sp_stats.t.sf(np.abs(t_stats), df=n - p)
        except np.linalg.LinAlgError:
            p_values = np.array([1.0, 1.0, 1.0])
    else:
        p_values = np.array([1.0, 1.0, 1.0])

    significant = []
    if p_values[0] < P_VALUE_THRESHOLD:
        significant.append('dipole')
    if p_values[1] < P_VALUE_THRESHOLD:
        significant.append('quadrupole')

    return {
        'A_dipole': A_dipole,
        'B_quadrupole': B_quadrupole,
        'intercept': intercept,
        'r_squared': r_squared,
        'p_dipole': p_values[0],
        'p_quadrupole': p_values[1],
        'rms': rms,
        'significant_terms': significant
    }


def generate_decoy_axes(n_decoys=N_DECOY_AXES, seed=DECOY_AXES_SEED):
    """
    Generate random sky axes uniform on the sphere.

    Parameters
    ----------
    n_decoys : int
        Number of decoy axes
    seed : int
        Random seed for reproducibility

    Returns
    -------
    numpy.ndarray
        Decoy unit vectors, shape (n_decoys, 3)
    """
    rng = np.random.RandomState(seed)

    ra = rng.uniform(0, 360, n_decoys)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n_decoys)))

    vectors = radec_to_unit_vector(ra, dec)

    return vectors, ra, dec


def run_single_axis_analysis(df, axis_vector):
    """
    Run complete dipole + quadrupole analysis for a single axis.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: ra, dec, delta_z
    axis_vector : numpy.ndarray
        Unit vector defining the test axis

    Returns
    -------
    dict
        Complete analysis results for this axis
    """
    cos_theta = compute_cos_theta(df['ra'].values, df['dec'].values, axis_vector)

    dipole_result = fit_dipole(cos_theta, df['delta_z'].values)
    dq_result = fit_dipole_quadrupole(cos_theta, df['delta_z'].values)

    return {
        'dipole_r': dipole_result['r'],
        'dipole_p': dipole_result['p_value'],
        'dipole_amplitude': dipole_result['A_dipole'],
        'dipole_rms': dipole_result['rms'],
        'dq_r_squared': dq_result['r_squared'],
        'dq_p_dipole': dq_result['p_dipole'],
        'dq_p_quadrupole': dq_result['p_quadrupole'],
        'dq_A_dipole': dq_result['A_dipole'],
        'dq_B_quadrupole': dq_result['B_quadrupole'],
        'dq_rms': dq_result['rms'],
        'significant_terms': dq_result['significant_terms']
    }


def run_blind_test(df, dataset_name="unknown"):
    """
    Run the complete blind test for one dataset:
    1. Compute residuals
    2. Analyse true axis (CMB dipole velocity direction)
    3. Analyse N_DECOY decoy axes
    4. Rank true axis among all axes
    5. Apply PASS/FAIL criteria

    Parameters
    ----------
    df : pd.DataFrame
        Filtered dataset with columns: ra, dec, z_obs
    dataset_name : str
        Name for reporting

    Returns
    -------
    dict
        Complete test results including verdict
    """
    df = compute_residuals(df)

    true_result = run_single_axis_analysis(df, V_HAT)

    decoy_vectors, decoy_ra, decoy_dec = generate_decoy_axes()
    decoy_results = []
    for i in range(len(decoy_vectors)):
        dr = run_single_axis_analysis(df, decoy_vectors[i])
        dr['axis_id'] = f'decoy_{i+1}'
        dr['axis_ra'] = decoy_ra[i]
        dr['axis_dec'] = decoy_dec[i]
        decoy_results.append(dr)

    all_abs_r = [abs(true_result['dipole_r'])]
    for dr in decoy_results:
        all_abs_r.append(abs(dr['dipole_r']))

    sorted_indices = np.argsort(all_abs_r)[::-1]
    true_rank = int(np.where(sorted_indices == 0)[0][0]) + 1

    passes_r = abs(true_result['dipole_r']) >= CORRELATION_THRESHOLD
    passes_p = true_result['dipole_p'] < P_VALUE_THRESHOLD
    passes_rank = true_rank == 1
    passes_all = passes_r and passes_p and passes_rank

    return {
        'dataset_name': dataset_name,
        'n_objects': len(df),
        'z_obs_mean': df['z_obs'].mean(),
        'z_model_range': [df['z_model'].min(), df['z_model'].max()],

        'true_r': true_result['dipole_r'],
        'true_p': true_result['dipole_p'],
        'true_amplitude': true_result['dipole_amplitude'],
        'true_rms': true_result['dipole_rms'],
        'true_dq_significant': true_result['significant_terms'],

        'decoy_r_values': [dr['dipole_r'] for dr in decoy_results],
        'decoy_p_values': [dr['dipole_p'] for dr in decoy_results],

        'true_rank': true_rank,
        'total_axes': 1 + len(decoy_results),

        'passes_r_threshold': passes_r,
        'passes_p_threshold': passes_p,
        'passes_rank': passes_rank,

        'VERDICT': 'PASS' if passes_all else 'FAIL'
    }
