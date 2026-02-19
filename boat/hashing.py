"""
BOAT Hashing Module
SHA-256 hashing for datasets, outputs, and manifest verification.
Ensures immutability of all test artifacts.
"""

import hashlib
import json
import os


def hash_file(filepath):
    """
    Compute SHA-256 hash of a file.

    Parameters
    ----------
    filepath : str
        Path to file

    Returns
    -------
    str
        Hex digest of SHA-256 hash
    """
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def hash_dataframe(df):
    """
    Compute SHA-256 hash of a pandas DataFrame.
    Converts to CSV bytes for deterministic hashing.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to hash

    Returns
    -------
    str
        Hex digest of SHA-256 hash
    """
    csv_bytes = df.to_csv(index=False).encode('utf-8')
    return hashlib.sha256(csv_bytes).hexdigest()


def hash_dict(d):
    """
    Compute SHA-256 hash of a dictionary.
    Serialises to sorted JSON for deterministic hashing.

    Parameters
    ----------
    d : dict
        Dictionary to hash (must be JSON-serialisable)

    Returns
    -------
    str
        Hex digest of SHA-256 hash
    """
    json_bytes = json.dumps(d, sort_keys=True, default=str).encode('utf-8')
    return hashlib.sha256(json_bytes).hexdigest()


def save_results_json(results, filepath):
    """
    Save results dictionary to JSON file.

    Parameters
    ----------
    results : dict
        Results to save
    filepath : str
        Output path

    Returns
    -------
    str
        SHA-256 hash of the saved file
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    return hash_file(filepath)


def generate_artifact_manifest(artifacts):
    """
    Generate a manifest of all artifact files with their hashes.

    Parameters
    ----------
    artifacts : dict
        {filename: filepath} mapping

    Returns
    -------
    dict
        {filename: sha256_hash} mapping
    """
    return {name: hash_file(path) for name, path in artifacts.items()
            if os.path.exists(path)}
