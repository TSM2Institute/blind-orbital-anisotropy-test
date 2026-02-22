"""
DESI DR1 data download via NOIRLab TAP query service.
Downloads BGS and LRG spectroscopic catalogues in redshift chunks.
"""

import urllib.request
import urllib.parse
import pandas as pd
import io
import os
import time

TAP_URL = "https://datalab.noirlab.edu/tap/sync"
OUTPUT_DIR = "attached_assets"


def tap_query(query, timeout=120):
    """Execute a TAP sync query and return CSV as string."""
    params = urllib.parse.urlencode({
        'QUERY': query,
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'FORMAT': 'csv'
    })
    url = f"{TAP_URL}?{params}"

    req = urllib.request.Request(url)
    resp = urllib.request.urlopen(req, timeout=timeout)
    return resp.read().decode()


def download_survey(survey_label, program, z_min, z_max, z_step, output_filename):
    """
    Download a DESI survey in redshift chunks.

    Args:
        survey_label: "BGS" or "LRG" (for logging)
        program: "bright" (BGS) or "dark" (LRG)
        z_min: minimum redshift
        z_max: maximum redshift
        z_step: chunk size
        output_filename: final CSV filename in OUTPUT_DIR
    """
    output_path = os.path.join(OUTPUT_DIR, output_filename)

    if os.path.exists(output_path):
        df_existing = pd.read_csv(output_path)
        print(f"{survey_label}: Already exists ({len(df_existing):,} rows). Skipping.")
        return output_path

    print(f"\n{'='*60}")
    print(f"Downloading DESI DR1 {survey_label}")
    print(f"Program: {program}, z range: {z_min}–{z_max}, step: {z_step}")
    print(f"{'='*60}")

    all_dfs = []
    z_lo = z_min
    chunk_num = 0

    while z_lo < z_max:
        z_hi = min(z_lo + z_step, z_max)
        chunk_num += 1

        query = f"""
        SELECT p.ra, p.dec, z.z
        FROM desi_dr1.zpix AS z
        JOIN desi_dr1.photometry AS p ON z.targetid = p.targetid
        WHERE z.zwarn = 0
            AND z.spectype = 'GALAXY'
            AND z.survey = 'main'
            AND z.program = '{program}'
            AND z.zcat_primary = 'True'
            AND z.z >= {z_lo}
            AND z.z < {z_hi}
        """

        print(f"  Chunk {chunk_num}: z = {z_lo:.3f}–{z_hi:.3f} ... ", end="", flush=True)

        t0 = time.time()
        try:
            csv_data = tap_query(query, timeout=300)
            df_chunk = pd.read_csv(io.StringIO(csv_data))
            elapsed = time.time() - t0
            print(f"{len(df_chunk):,} rows ({elapsed:.1f}s)")

            if len(df_chunk) > 0:
                all_dfs.append(df_chunk)
        except Exception as e:
            elapsed = time.time() - t0
            print(f"ERROR ({elapsed:.1f}s): {e}")
            print(f"  Retrying chunk {chunk_num} in 10 seconds...")
            time.sleep(10)
            try:
                csv_data = tap_query(query, timeout=300)
                df_chunk = pd.read_csv(io.StringIO(csv_data))
                print(f"  Retry successful: {len(df_chunk):,} rows")
                if len(df_chunk) > 0:
                    all_dfs.append(df_chunk)
            except Exception as e2:
                print(f"  Retry also failed: {e2}")
                print(f"  SKIPPING chunk {chunk_num}. Manual retry may be needed.")

        z_lo = z_hi

        time.sleep(1)

    if not all_dfs:
        print(f"ERROR: No data downloaded for {survey_label}")
        return None

    df_all = pd.concat(all_dfs, ignore_index=True)

    df_all.columns = ['ra', 'dec', 'z']

    n_before = len(df_all)
    df_all = df_all.drop_duplicates()
    n_dupes = n_before - len(df_all)

    df_all.to_csv(output_path, index=False)

    print(f"\n{survey_label} download complete:")
    print(f"  Total galaxies: {len(df_all):,}")
    print(f"  Duplicates removed: {n_dupes:,}")
    print(f"  z range: {df_all['z'].min():.4f} – {df_all['z'].max():.4f}")
    print(f"  Median z: {df_all['z'].median():.4f}")
    print(f"  Saved to: {output_path}")

    return output_path


if __name__ == "__main__":
    download_survey(
        survey_label="BGS",
        program="bright",
        z_min=0.001,
        z_max=0.600,
        z_step=0.05,
        output_filename="DESI_DR1_BGS.csv"
    )

    download_survey(
        survey_label="LRG",
        program="dark",
        z_min=0.300,
        z_max=1.100,
        z_step=0.05,
        output_filename="DESI_DR1_LRG.csv"
    )

    print("\n\nDone. Both catalogues saved to attached_assets/")
    print("Next: update datasets.py DESI loaders, then run MC calibration.")
