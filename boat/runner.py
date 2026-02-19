"""
BOAT Runner — Main Test Orchestrator
Executes the complete Blind Orbital Anisotropy Test.
Produces final PASS/FAIL verdict and all immutable artifacts.
"""

import os
import json
import time
import numpy as np
import pandas as pd
from datetime import datetime, timezone

from boat.manifest import (
    N_DATASETS, PASS_REQUIRED, CORRELATION_THRESHOLD,
    P_VALUE_THRESHOLD, V_OBS, CST_PERIOD_GYR,
    COG_X_RA_DEG, COG_X_DEC_DEG, CMB_DIPOLE_RA_DEG, CMB_DIPOLE_DEC_DEG,
    R_OBS_MLY, BETA_OBS, DATASET_SELECTION_SEED, DECOY_AXES_SEED,
    N_DECOY_AXES, Z_MIN, Z_MAX, GALACTIC_PLANE_MASK_DEG, MIN_SAMPLE_SIZE
)
from boat.datasets import load_dataset, apply_quality_cuts, select_datasets
from boat.analysis import run_blind_test
from boat.hashing import hash_dataframe, hash_dict, save_results_json


RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')


def run_full_test(mode='synthetic'):
    """
    Execute the complete BOAT test.

    Parameters
    ----------
    mode : str
        'synthetic' for validation run, 'production' for real data

    Returns
    -------
    dict
        Complete test results with overall verdict
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("  BLIND ORBITAL ANISOTROPY TEST (BOAT)")
    print("  TSM2 Institute for Cosmology Ltd")
    print(f"  Execution time: {datetime.now(timezone.utc).isoformat()}")
    print(f"  Mode: {mode}")
    print("=" * 70)

    manifest_record = {
        'v_obs_km_s': V_OBS,
        'cst_period_gyr': CST_PERIOD_GYR,
        'cog_x_ra_deg': COG_X_RA_DEG,
        'cog_x_dec_deg': COG_X_DEC_DEG,
        'cmb_dipole_ra_deg': CMB_DIPOLE_RA_DEG,
        'cmb_dipole_dec_deg': CMB_DIPOLE_DEC_DEG,
        'r_obs_mly': R_OBS_MLY,
        'beta_obs': BETA_OBS,
        'correlation_threshold': CORRELATION_THRESHOLD,
        'p_value_threshold': P_VALUE_THRESHOLD,
        'n_datasets': N_DATASETS,
        'pass_required': PASS_REQUIRED,
        'n_decoy_axes': N_DECOY_AXES,
        'z_min': Z_MIN,
        'z_max': Z_MAX,
        'galactic_mask_deg': GALACTIC_PLANE_MASK_DEG,
        'min_sample_size': MIN_SAMPLE_SIZE,
        'dataset_selection_seed': DATASET_SELECTION_SEED,
        'decoy_axes_seed': DECOY_AXES_SEED,
    }

    print("\n--- Dataset Selection ---")
    datasets = select_datasets()
    print(f"  Selected {len(datasets)} datasets: {datasets}")

    all_results = []
    dataset_hashes = {}

    for i, ds_name in enumerate(datasets):
        print(f"\n{'='*70}")
        print(f"  DATASET {i+1}/{len(datasets)}: {ds_name}")
        print(f"{'='*70}")

        if mode == 'synthetic':
            seed = 42 + i
        else:
            seed = DATASET_SELECTION_SEED + i
        df_raw = load_dataset(ds_name, seed=seed)
        print(f"  Raw objects: {len(df_raw)}")

        df_filtered, filter_stats = apply_quality_cuts(df_raw, ds_name)
        print(f"  After quality cuts: {filter_stats['n_final']}")

        if not filter_stats['passes_min_size']:
            print(f"  WARNING: Below minimum sample size ({MIN_SAMPLE_SIZE})")
            result = {
                'dataset_name': ds_name,
                'n_objects': filter_stats['n_final'],
                'VERDICT': 'FAIL',
                'reason': 'Below minimum sample size',
                'filter_stats': filter_stats
            }
        else:
            dataset_hashes[ds_name] = hash_dataframe(df_filtered)
            print(f"  Dataset hash: {dataset_hashes[ds_name][:16]}...")

            result = run_blind_test(df_filtered, ds_name)
            result['filter_stats'] = filter_stats
            result['dataset_hash'] = dataset_hashes[ds_name]

        print(f"\n  --- Results for {ds_name} ---")
        if 'true_r' in result:
            print(f"  True axis |r|: {abs(result['true_r']):.6f} "
                  f"(threshold: {CORRELATION_THRESHOLD})")
            print(f"  True axis p:   {result['true_p']:.4e} "
                  f"(threshold: {P_VALUE_THRESHOLD})")
            print(f"  True axis rank: {result['true_rank']} of "
                  f"{result['total_axes']}")

            if 'decoy_r_values' in result:
                decoy_abs_r = [abs(r) for r in result['decoy_r_values']]
                print(f"  Decoy |r| range: [{min(decoy_abs_r):.6f}, "
                      f"{max(decoy_abs_r):.6f}]")

            print(f"  Passes |r|: {result['passes_r_threshold']}")
            print(f"  Passes p:   {result['passes_p_threshold']}")
            print(f"  Passes rank: {result['passes_rank']}")

        print(f"  >>> VERDICT: {result['VERDICT']} <<<")

        all_results.append(result)

    n_pass = sum(1 for r in all_results if r['VERDICT'] == 'PASS')
    n_fail = sum(1 for r in all_results if r['VERDICT'] == 'FAIL')
    overall_verdict = 'PASS' if n_pass >= PASS_REQUIRED else 'FAIL'

    output = {
        'test_name': 'Blind Orbital Anisotropy Test (BOAT)',
        'institution': 'TSM2 Institute for Cosmology Ltd',
        'execution_time': datetime.now(timezone.utc).isoformat(),
        'mode': mode,
        'manifest_constants': manifest_record,
        'datasets_tested': len(all_results),
        'n_pass': n_pass,
        'n_fail': n_fail,
        'pass_required': PASS_REQUIRED,
        'overall_verdict': overall_verdict,
        'per_dataset_results': all_results,
        'dataset_hashes': dataset_hashes,
    }

    results_path = os.path.join(RESULTS_DIR, f'boat_results_{mode}.json')
    results_hash = save_results_json(output, results_path)
    output['results_file_hash'] = results_hash

    print("\n" + "=" * 70)
    print("  FINAL RESULTS")
    print("=" * 70)
    print(f"  Datasets tested: {len(all_results)}")
    print(f"  Passed: {n_pass}")
    print(f"  Failed: {n_fail}")
    print(f"  Required to pass: {PASS_REQUIRED}")
    print()
    print(f"  ╔══════════════════════════════════════╗")
    print(f"  ║  OVERALL VERDICT:  {overall_verdict:^17s} ║")
    print(f"  ╚══════════════════════════════════════╝")
    print()
    print(f"  Results saved: {results_path}")
    print(f"  Results hash:  {results_hash[:16]}...")
    print("=" * 70)

    return output


if __name__ == '__main__':
    run_full_test(mode='synthetic')
