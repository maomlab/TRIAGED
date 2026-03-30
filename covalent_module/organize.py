'''
Boltz Prediction Reorganization Script

Standalone script to reorganize Boltz-2 predictions after docking jobs complete.
Uses the same JSON configuration format as main.py.

Workflow:
---------
1. Parse JSON config (same format as main.py)
2. Update errored.csv with compounds missing .cif predictions
3. Update predictions.csv with confidence metrics
4. Reorganize CIF + YAML files into structured output directory
5. Remove temporary PKL files (covalent mode only)
6. Clean up run_cache directory

Usage:
------
    python reorg_predictions.py --json_file config.json

JSON Configuration Format:
--------------------------
{
    "RECORD_PATH": "/path/to/records",
    "PDB": "/path/to/protein.pdb",
    "RUN_CACHE": "/path/to/temp/cache",
    "BOLTZ_CACHE": "/path/to/boltz/weights",
    "OUTPUT": "/path/to/output",
    "COVALENT": true,
    "VERBOSE": true
}

Authors: Manasa Yadavalli
'''

import os
import json
import glob
import shutil
import argparse
import pandas as pd

from BoltzCov.update_ligands import update_predictions


def read_json_args(json_file):
    '''
    Reads the JSON input file with arguments required for prediction reorganization.

    :param json_file: Path to JSON file.
    :return: Tuple of parsed configuration values.
    '''
    with open(json_file, 'r') as jf:
        run_arguments = json.load(jf)

    record_path = run_arguments.get("RECORD_PATH")
    record_path = os.path.expandvars(record_path) if record_path else None

    pdb = run_arguments.get("PDB")
    pdb = os.path.expandvars(pdb) if pdb else None

    run_cache = run_arguments.get("RUN_CACHE")
    run_cache = os.path.expandvars(run_cache) if run_cache else None

    boltz_cache = run_arguments.get("BOLTZ_CACHE")
    boltz_cache = os.path.expandvars(boltz_cache) if boltz_cache else None

    output_dir = run_arguments.get("OUTPUT")
    output_dir = os.path.expandvars(output_dir) if output_dir else None

    COVALENT = run_arguments.get("COVALENT", False)
    VERBOSE = run_arguments.get("VERBOSE", False)

    return record_path, pdb, run_cache, boltz_cache, output_dir, COVALENT, VERBOSE


def main(args):
    '''
    Reorganizes Boltz-2 predictions using the same JSON config as main.py.

    Steps:
    1. Parse config and validate required arguments
    2. Update errored.csv and predictions.csv from run_cache outputs
    3. Reorganize CIF/YAML files into output_dir
    4. Remove temporary PKL files (covalent only)
    5. Delete run_cache directory
    '''
    record_path, pdb, run_cache, boltz_cache, output_dir, COVALENT, VERBOSE = read_json_args(args.json_file)

    # Validate required arguments
    required = {
        'RECORD_PATH': record_path,
        'PDB': pdb,
        'RUN_CACHE': run_cache,
        'BOLTZ_CACHE': boltz_cache,
        'OUTPUT': output_dir,
    }

    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")

    protein_name = os.path.splitext(os.path.basename(pdb))[0]
    run_cache_prot = os.path.join(run_cache, protein_name)

    if not os.path.exists(run_cache_prot):
        raise FileNotFoundError(
            f"Run cache directory not found: {run_cache_prot}\n"
            f"Ensure docking jobs completed and run_cache/protein_name exists."
        )

    # Determine correct predictions CSV based on covalent mode
    if COVALENT:
        pred_csv = os.path.join(record_path, 'predictions.csv')
    else:
        pred_csv = os.path.join(record_path, 'noncov_predictions.csv')

    # Step 1: Update errored.csv and predictions CSV with results from run_cache
    print(f"1. Checking predictions and updating records in {record_path}...")

    if not COVALENT:
        # make temp ligand_csv 
        lig_dirs = os.listdir(run_cache_prot)
        lig_dirs = [d for d in os.listdir(run_cache_prot) if d[0].isdigit()]
        substance_ids = [lig.split('_')[-1] for lig in lig_dirs if os.path.isdir(os.path.join(run_cache_prot, lig))]
        ligand_noncov_df = pd.DataFrame({'pkl_id': substance_ids})
        ligand_noncov_df.to_csv(f'{run_cache_prot}/noncov_ligands.csv', index=False)

    update_predictions.check_pred(
        run_cache_prot=run_cache_prot,
        record_path=record_path,
        VERBOSE=VERBOSE,
        protein_name=protein_name,
        COVALENT=COVALENT
    )

    # Step 2: Load updated predictions and filter to current protein
    if not os.path.exists(pred_csv):
        raise FileNotFoundError(
            f"Predictions CSV not found after check_pred: {pred_csv}"
        )

    csv_files = glob.glob(os.path.join(run_cache_prot, '*.csv'))
    if csv_files: ligand_cache_csv = csv_files[0] # temp csv that was generated before boltz run to make yamls 
    ligand_csv_df = pd.read_csv(ligand_cache_csv)

    if VERBOSE:
        print(f"   Found {len(ligand_csv_df)} predictions for {protein_name}.")

    # Step 3: Reorganize CIF + YAML files into output directory
    # Skip compounds whose CIF already exists in temp output_dir 
    output_temp = os.path.join(run_cache_prot, 'temp', protein_name) # should be unique per run

    print(f"2. Reorganizing {len(ligand_csv_df)} predictions into {output_temp}...")
    update_predictions.reorg_preds(run_cache_prot, ligand_csv_df, output_temp, VERBOSE)
    
    # Step 4: Remove temporary PKL files (covalent mode only)
    if COVALENT:
        print("3. Removing temporary PKL files (covalent mode)...")
        update_predictions.remove_pkls(boltz_cache, run_cache_prot, VERBOSE)

    # Step 5: Reorg into replicate dirs
    if protein_name == '3F75':
        protein= 'tgcpl' 
    elif protein_name == '5MAJ':
        protein = 'hscpl'

    update_predictions.reorg_reps(output_temp, output_dir, protein)

    # Step 6: Clean up run_cache with permission
    run_cache_parent = os.path.dirname(run_cache_prot)
    print(f"4. Run cache directory: {run_cache_parent}")
    confirm = input("Delete run cache directory? This cannot be undone. [y/n]: ").strip().lower()
    if confirm == 'y':
        shutil.rmtree(run_cache_parent)
        print(f"   Deleted {run_cache_parent}.")
    else:
        print("   Skipping deletion. Run cache preserved.")

    print("Done. Prediction reorganization complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reorganize Boltz-2 predictions using a JSON config file."
    )
    parser.add_argument(
        "--json_file",
        required=True,
        help="Path to JSON configuration file (same format as main.py)."
    )
    args = parser.parse_args()
    main(args)