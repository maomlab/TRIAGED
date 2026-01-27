import os
import pandas as pd
import sys
import json 
import time 
import argparse
from BoltzCov.update_ligands import pull_cdd_ligs

def read_json_args(json_file):
    '''
    Reads the json input file with arguments required to dock covalently with Boltz-2.
    
    :param json_file: Path to json file.
    '''

    with open(json_file, 'r') as jf:
        run_arguments = json.load(jf)

        record_path = run_arguments.get("RECORD_PATH") 
        record_path = os.path.expandvars(record_path) if record_path else None  # can be None for first time run since ligands are pulled from cdd but path must be given 

        cdd_api_key = run_arguments.get("CDD_API_KEY", None)
        vault_id = run_arguments.get("VAULT_ID", None) # number of the vault in CDD 
        readout_query = run_arguments.get("READOUT_KEY", None) # is a dict 
        mol_query = run_arguments.get("MOL_QUERY", None) # is a dict 
        
        pdb = run_arguments.get("PDB")
        pdb = os.path.expandvars(pdb) if pdb else None 
        res_idx = run_arguments.get("RES_IDX", None)
        ligand_chain = run_arguments.get("LIGAND_CHAIN", None)
        msa_path = run_arguments.get("MSA_PATH")
        msa_path = os.path.expandvars(msa_path) if msa_path else None # dont need to provide this  
        boltz_cache = run_arguments.get("BOLTZ_CACHE", None)
        boltz_cache = os.path.expandvars(boltz_cache) if boltz_cache else None 
        slurm_template = run_arguments.get("SLURM_TEMPLATE", None)
        slurm_template = os.path.expandvars(slurm_template) if slurm_template else None 
        
        VERBOSE = run_arguments.get("VERBOSE", True)
        output_dir = run_arguments.get("OUTPUT", None)
        output_dir = os.path.expandvars(output_dir) if output_dir else None 
    
    return (
    record_path,
    (cdd_api_key, vault_id, readout_query, mol_query),
    (pdb, res_idx, ligand_chain, msa_path, boltz_cache, slurm_template),
    (VERBOSE, output_dir)
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
    Future features:
    6) takes either cdd api or ligand csv 
    7) pulls replicate readouts as well 
    8) Produces analysis plots 
    9) Computes interaction fingerprints and interactions 
    10) Produces visualizations
    '''
    # loading input arguments from user 
    (record_path,
    (cdd_api_key, vault_id, readout_query, mol_query, readout_id),
    (pdb, res_idx, ligand_chain, msa_path, boltz_cache, slurm_template),
    (VERBOSE, output_dir)
    ) = read_json_args(args.json_file)

    if any(x is None for x in (pdb, res_idx, ligand_chain, boltz_cache, slurm_template,
    cdd_api_key, vault_id, readout_query, mol_query, readout_id, output_dir)):
        raise ValueError("Please make sure all required arguments are given in the input JSON.")

    # pull ligand information and experiment readouts 
    print("1. Pulling ligands from CDD vault using the following queries:" \
        f"{readout_query}" \
        f"{mol_query}")
            
    readouts, molecules = pull_cdd_ligs.cdd_query(API_KEY=cdd_api_key, 
                                    VAULT_ID=7171, 
                                    readout_query=readout_query, 
                                    mol_query=mol_query)
    
    print("2. Updating compound records.")
    new_readouts_df, new_metadata_df = pull_cdd_ligs.get_ic50s(readouts, molecules)
    if record_path is None: 
        if VERBOSE: print("Writing new records (metadata.csv, and experimental_readouts.csv) since None path provided by User.")
        # make record dir in output if dir dne 
        record_dir = os.path.join(output_dir, 'records')
        os.makedirs(record_dir, exist_ok=True)

        pd.DataFrame(new_readouts_df).to_csv(os.path.join(record_dir,"experiment_readouts.csv"), index=False)
        pd.DataFrame(new_metadata_df).to_csv(os.path.join(record_dir, "metadata.csv"), index=False)

    if record_path: 
        os.makedirs(record_path, exist_ok=True)
        old_meta = os.path.join(record_path, 'metadata.csv') 
        old_exp = os.path.join(record_path, 'experiment_readouts.csv')
        if not os.path.exists(old_meta) or not os.path.exists(old_exp):
            if VERBOSE: print(f"metadata.csv or/and experiment_readouts.csv were not found in {record_path}. Writing new records.") # first time use or to skirt unintential overwriting 
            print("You have 5 seconds to terminate and cancel overwrite to possible exisiting records.")
            time.sleep(5)
            pd.DataFrame(new_readouts_df).to_csv(old_exp, index=False) # writing new records to record path provided by user
            pd.DataFrame(new_metadata_df).to_csv(old_meta, index=False)
            if VERBOSE: print(f"Fresh metadata.csv and experiment_readouts.csv written in {record_path}.")

        if os.path.exists(old_meta) and os.path.exists(old_exp):
            if VERBOSE: print(f"Updating the provided metadata.csv and experiment_readouts.csv in {record_path}")
            old_meta_df = pd.read_csv(old_meta)
            old_exp_df = pd.read_csv(old_exp)
            # only needs to be updated in cases where we have old record files existing under same names 
            updated_metadata, updated_readouts = pull_cdd_ligs.update_local_data(old_meta_df, old_exp_df, new_metadata_df, new_readouts_df)
            # overwriting existing record files to update
            pd.DataFrame(updated_metadata).to_csv(old_meta, index=False)
            pd.DataFrame(updated_readouts).to_csv(old_exp, index=False)
    
    if VERBOSE: print("-SUCCESS- 2. Record files updated with CDD-Vault information.")

    # CHECK IF PREDICTION CSV EXISTS BEFORE UPDATING IT 
    print("3. Checking exisiting predictions and performing Docking with Boltz-2 Covalent.")
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_file", required=True)
    args = parser.parse_args()

    main(args)