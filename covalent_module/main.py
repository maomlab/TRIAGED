'''
Boltz Covalent Docking Pipeline - Main Workflow Script
Env: cdd_pkl
This script orchestrates the complete workflow for covalent and noncovalent docking using Boltz-2, 
integrating with CDD Vault for compound management and experimental data tracking.

Workflow Overview:
-----------------
1. Query CDD Vault for ligands and experimental readouts
2. Update or create local compound records (metadata.csv, experiment_readouts.csv)
3. Identify compounds that need docking (exclude previously attempted)
4. Submit Boltz-2 covalent docking jobs via SLURM
5. Monitor job completion

Input Requirements:
------------------
- JSON configuration file with parameters

Output Structure:
----------------
<OUTPUT_DIR>/
├── <replicate_1>/
│   ├── PROT_LIG_model_0.cif
│   └── hparams.yaml
├── <replicate_2>/
│   └── ...
└── records/
    ├── metadata.csv               # Compound metadata from CDD
    ├── experiment_readouts.csv    # Experimental IC50/readout data
    ├── predictions.csv            # Boltz prediction metrics (compiled)
    └── errored.csv               # Failed compounds (to skip on retry)

JSON Configuration Format:
-------------------------
{
    "RECORD_PATH": "/path/to/records",
    "CDD_API_KEY": "your_api_key",
    "VAULT_ID": 12345,
    "READOUT_QUERY": {"protocol_ids": [123], "runs": [456]},
    "MOL_QUERY": {"molecule_ids": [789]},
    "PDB": "/path/to/protein.pdb",
    "RES_IDX": 145,
    "LIGAND_CHAIN": "X",
    "MSA_PATH": "/path/to/msa/dir",
    "RUN_CACHE": "/path/to/temp/cache",
    "BOLTZ_CACHE": "/path/to/boltz/weights",
    "SLURM_TEMPLATE": "/path/to/slurm_template.sh",
    "VERBOSE": true,
    "OUTPUT": "/path/to/output"
}

Usage:
------
    python main.py --json_file config.json

Dependencies:
------------
- BoltzCov package (update_ligands, run_boltz modules)
- pandas, json, shutil, time, argparse, os
- CDD Vault API access
- SLURM cluster environment
- Boltz-2 model weights

Notes:
-----
- run_cache is deleted and recreated on each run (5 second warning)
- Boltz weights are downloaded automatically to boltz_cache if missing
- predictions.csv is incrementally updated (preserves previous results)
- Failed docking attempts are logged to errored.csv to avoid retries

Authors: Manasa Yadavalli
'''
import os
import pandas as pd
import json 
import shutil
import time 
import argparse
from BoltzCov.update_ligands import pull_cdd_ligs, update_predictions
from BoltzCov.run_boltz import pull_boltz2_weights

# use ccd_pkl
def read_json_args(json_file):
    '''
    Reads the json input file with arguments required to dock covalently with Boltz-2.
    
    :param json_file: Path to json file.
    '''

    with open(json_file, 'r') as jf:
        run_arguments = json.load(jf)

        record_path = run_arguments.get("RECORD_PATH") # can have csvs named: metadata.csv, experiment_readout.csv, predictions.csv
        record_path = os.path.expandvars(record_path) if record_path else None  # can be None for first time run since ligands are pulled from cdd

        cdd_api_key = run_arguments.get("CDD_API_KEY", None)
        vault_id = run_arguments.get("VAULT_ID", None) # number of the vault in CDD 
        readout_query = run_arguments.get("READOUT_QUERY", None) # is a dict 
        mol_query = run_arguments.get("MOL_QUERY", None) # is a dict 
        syn_include = run_arguments.get("SYN_INCLUDE", None)
        syn_exclude = run_arguments.get("SYN_EXCLUDE", None)
        
        pdb = run_arguments.get("PDB")
        pdb = os.path.expandvars(pdb) if pdb else None 
        res_idx = run_arguments.get("RES_IDX", None)
        ligand_chain = run_arguments.get("LIGAND_CHAIN", None)
        msa_path = run_arguments.get("MSA_PATH")
        msa_path = os.path.expandvars(msa_path) if msa_path else None # dont need to provide this  
        run_cache = run_arguments.get("RUN_CACHE", None)
        run_cache = os.path.expandvars(run_cache) if run_cache else None  # need this; just easier to have a single location to put temp stuff in and delete it completely after the run
        boltz_cache = run_arguments.get("BOLTZ_CACHE", None)
        boltz_cache = os.path.expandvars(boltz_cache) if boltz_cache else None  # this is needed for mols/weights downloading 
        slurm_template = run_arguments.get("SLURM_TEMPLATE", None)
        slurm_template = os.path.expandvars(slurm_template) if slurm_template else None 
        
        VERBOSE = run_arguments.get("VERBOSE", False)
        output_dir = run_arguments.get("OUTPUT", None) # include protein name if you want it to be stored in a seperate protein directory 
        output_dir = os.path.expandvars(output_dir) if output_dir else None 
        COVALENT = run_arguments.get("COVALENT", False)
    return (
    record_path,
    cdd_api_key, vault_id, readout_query, mol_query, syn_include, syn_exclude,
    pdb, res_idx, ligand_chain, msa_path, boltz_cache, run_cache, slurm_template,
    VERBOSE, output_dir, COVALENT
    )

def main(args):
    '''
    Takes in JSON with all required variables.
    Top prority: 
    1) Pulls ligands from CDD-Vault based on queries
    2) Updates compound records
    3) Docks ligands that we did not dock previously 
    4) Reorginizes outputs 
    5) Updates relavant results in compound records
    '''
    # loading input arguments from user 
    (record_path,
    cdd_api_key, vault_id, readout_query, mol_query, syn_include, syn_exclude,
    pdb, res_idx, ligand_chain, msa_path, boltz_cache, run_cache, slurm_template,
    VERBOSE, output_dir, COVALENT) = read_json_args(args.json_file)

    if COVALENT:
        missing = []
        cov_params = {
            'pdb': pdb, 
            'res_idx': res_idx, 
            'ligand_chain': ligand_chain, 
            'slurm_template': slurm_template,
            'cdd_api_key': cdd_api_key, 
            'vault_id': vault_id, 
            'readout_query': readout_query, 
            'mol_query': mol_query, 
            'boltz_cache': boltz_cache, 
            'run_cache':run_cache,
            'output_dir': output_dir
        }

        for name, value in cov_params.items():
            if value is None:
                missing.append(name)
    else:
        missing = []
        params = {
            'pdb': pdb, 
            'slurm_template': slurm_template,
            'cdd_api_key': cdd_api_key, 
            'vault_id': vault_id, 
            'readout_query': readout_query, 
            'mol_query': mol_query, 
            'boltz_cache': boltz_cache, 
            'run_cache':run_cache,
            'output_dir': output_dir
        }

        for name, value in params.items():
            if value is None:
                missing.append(name)

    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")

    if os.path.exists(run_cache):
        print("[WARNING] Run Cache exists and will be deleted.")
        print("Cancel in 5 seconds to prevent deletion. Ensure the cache directory is empty or choose a different directory.")
        time.sleep(5)
        shutil.rmtree(run_cache)
    os.makedirs(run_cache)

    if not os.path.exists(boltz_cache):
        os.makedirs(boltz_cache)
        from pathlib import Path
        print("[WARNING] Boltz weights do not exists.")
        print("Downloading weights and PDB ligand files...")
        pull_boltz2_weights.download_boltz2(Path(boltz_cache))

    # pull ligand information and experiment readouts 
    print("1. Pulling ligands from CDD vault using the following queries:\n"
      f"Readout query: {readout_query}\n"
      f"Molecule query: {mol_query}")
    
    readouts, molecules = pull_cdd_ligs.cdd_query(API_KEY=cdd_api_key, 
                                    VAULT_ID=vault_id,  # Use vault_id variable
                                    readout_query=readout_query, 
                                    mol_query=mol_query)
    
    print("2. Updating compound records.")
    old_meta = os.path.join(record_path, 'metadata.csv') 
    old_exp = os.path.join(record_path, 'experiment_readouts.csv')

    new_readouts_df, new_metadata_df = pull_cdd_ligs.get_ic50s(readouts, molecules, syn_include=syn_include, syn_exclude=syn_exclude)
    if record_path is None: 
        if VERBOSE: print("Writing new records (metadata.csv, and experimental_readouts.csv) in output directory since None path provided by User.")
        # make record dir in output if dir dne 
        record_dir = os.path.join(output_dir, 'records')
        os.makedirs(record_dir, exist_ok=True)

        pd.DataFrame(new_readouts_df).to_csv(os.path.join(record_dir,"experiment_readouts.csv"), index=False)
        pd.DataFrame(new_metadata_df).to_csv(os.path.join(record_dir, "metadata.csv"), index=False)

    if record_path: 
        os.makedirs(record_path, exist_ok=True)
        if not os.path.exists(old_meta) or not os.path.exists(old_exp):
            if VERBOSE: print(f"metadata.csv or/and experiment_readouts.csv were not found in {record_path}. Writing new records.") # first time use or to skirt unintential overwriting 
            print("You have 5 seconds to terminate and cancel overwrite to possible exisiting records.")
            time.sleep(5)
            pd.DataFrame(new_readouts_df).to_csv(old_exp, index=False) # writing new records to record path provided by user
            pd.DataFrame(new_metadata_df).to_csv(old_meta, index=False)
            if VERBOSE: print(f"Fresh metadata.csv and experiment_readouts.csv written in {record_path}.")

        elif os.path.exists(old_meta) and os.path.exists(old_exp):
            print("You have 5 seconds to terminate and cancel overwrite to possible exisiting records.")
            time.sleep(5)
            if VERBOSE: print(f"Updating the provided metadata.csv and experiment_readouts.csv in {record_path}")
            old_meta_df = pd.read_csv(old_meta)
            old_exp_df = pd.read_csv(old_exp)
            # only needs to be updated in cases where we have old record files existing under same names 
            updated_metadata, updated_readouts = pull_cdd_ligs.update_local_data(old_meta_df, old_exp_df, new_metadata_df, new_readouts_df)
            # overwriting existing record files to update
            if VERBOSE:
                new_meta_rows = len(updated_metadata) - len(old_meta_df)
                new_readout_rows = len(updated_readouts) - len(old_exp_df)
                print(f"Update will add ~{new_meta_rows} metadata and ~{new_readout_rows} readouts")
                confirm = input("Proceed? (y/n): ")
                if confirm.lower() != 'y':
                    print("Update cancelled")
                    raise ValueError('Execution Cancelled')
            
            pd.DataFrame(updated_metadata).to_csv(old_meta, index=False)
            pd.DataFrame(updated_readouts).to_csv(old_exp, index=False)
    
    if VERBOSE: print("-SUCCESS- Record files updated with CDD-Vault information.")
    
    # obtain list of ligands that need to be docked 
    protein_name = os.path.splitext(os.path.basename(pdb))[0]
    if COVALENT: 
        print("3. Checking existing predictions and performing Docking with Boltz-2 Covalent.")
        pred_rec = os.path.join(record_path, 'predictions.csv')
    else: 
        print("3. Checking existing predictions and performing Docking with Boltz-2 Non-Covalent.")
        pred_rec = os.path.join(record_path, 'noncov_predictions.csv') 
    metadata_df = pd.read_csv(old_meta)
    
    # exclude previous ligands that had a docking attempt 
    error_csv = os.path.join(record_path, 'errored.csv')
    if os.path.exists(pred_rec):
        pred_df = pd.read_csv(pred_rec)
        dock_compounds = update_predictions.fetch_new(pred_df, metadata_df, protein_name)
        if os.path.exists(error_csv):
            error_df = pd.read_csv(error_csv)
            dock_compounds = update_predictions.check_attempted(error_df, dock_compounds) 
        else: 
            if VERBOSE: print("No errored compounds.")

    else:
        if VERBOSE: print("predictions.csv was not found. Attempting to Dock all compounds.")
        dock_compounds = metadata_df[['substance_id', 'inchi_key', 'smiles']]
    
    if VERBOSE: 
        print(f"Docking {len(dock_compounds)} compounds. Cancel in 5 seconds to abort.")
        time.sleep(5)

    # submit jobs: run boltz    
    from BoltzCov.run_boltz import submit_job
    if COVALENT:
        final_status = submit_job.run_boltz_cov(prot_file=pdb, ligand_df=dock_compounds, boltz_cache=boltz_cache, run_cache=run_cache,
                        res_idx=res_idx, ligand_chain=ligand_chain, VERBOSE=VERBOSE, 
                        slurm_template=slurm_template, msa_path=msa_path)
    else:
        final_status = submit_job.run_boltz_noncov(prot_file=pdb, ligand_df=dock_compounds, boltz_cache=boltz_cache, 
                        run_cache=run_cache, VERBOSE=VERBOSE, slurm_template=slurm_template, msa_path=msa_path)
    # check when done 
    if final_status == "COMPLETED":
        print("Jobs completed!")
    else:
        print(f"Job failed with status: {final_status}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_file", required=True)
    args = parser.parse_args()

    main(args)