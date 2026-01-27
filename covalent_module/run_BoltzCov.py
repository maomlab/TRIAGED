import os
import pandas as pd
import json 
import shutil
import time 
import argparse
from BoltzCov.update_ligands import pull_cdd_ligs, update_predictions
from BoltzCov.run_boltz import submit_job

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
        readout_query = run_arguments.get("READOUT_KEY", None) # is a dict 
        mol_query = run_arguments.get("MOL_QUERY", None) # is a dict 
        
        pdb = run_arguments.get("PDB")
        pdb = os.path.expandvars(pdb) if pdb else None 
        res_idx = run_arguments.get("RES_IDX", None)
        ligand_chain = run_arguments.get("LIGAND_CHAIN", None)
        msa_path = run_arguments.get("MSA_PATH")
        msa_path = os.path.expandvars(msa_path) if msa_path else None # dont need to provide this  
        boltz_cache = run_arguments.get("BOLTZ_CACHE", None) 
        boltz_cache = os.path.expandvars(boltz_cache) if boltz_cache else None # need this; just easier to have a single location to put temp stuff in and delete it completely after the run
        slurm_template = run_arguments.get("SLURM_TEMPLATE", None)
        slurm_template = os.path.expandvars(slurm_template) if slurm_template else None 
        
        VERBOSE = run_arguments.get("VERBOSE", True)
        output_dir = run_arguments.get("OUTPUT", None)
        output_dir = os.path.expandvars(output_dir) if output_dir else None 
    
    return (
    record_path,
    cdd_api_key, vault_id, readout_query, mol_query,
    pdb, res_idx, ligand_chain, msa_path, boltz_cache, slurm_template,
    VERBOSE, output_dir
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
    cdd_api_key, vault_id, readout_query, mol_query, readout_id,
    pdb, res_idx, ligand_chain, msa_path, boltz_cache, slurm_template,
    VERBOSE, output_dir) = read_json_args(args.json_file)

    if any(x is None for x in (pdb, res_idx, ligand_chain, slurm_template,
    cdd_api_key, vault_id, readout_query, mol_query, boltz_cache, readout_id, output_dir)):
        raise ValueError("Please make sure all required arguments are given in the input JSON.")

    if os.path.exists(boltz_cache):
        print("[WARNING] Boltz Cache exists and will be deleting.")
        print("Cancel in 5 seconds to prevent deletion. Ensure the cache directory is empty.")
        time.sleep(5)
        shutil.rmtree(boltz_cache)
    os.makedirs(boltz_cache)

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

    new_readouts_df, new_metadata_df = pull_cdd_ligs.get_ic50s(readouts, molecules)
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
            if VERBOSE: print(f"Updating the provided metadata.csv and experiment_readouts.csv in {record_path}")
            old_meta_df = pd.read_csv(old_meta)
            old_exp_df = pd.read_csv(old_exp)
            # only needs to be updated in cases where we have old record files existing under same names 
            updated_metadata, updated_readouts = pull_cdd_ligs.update_local_data(old_meta_df, old_exp_df, new_metadata_df, new_readouts_df)
            # overwriting existing record files to update
            pd.DataFrame(updated_metadata).to_csv(old_meta, index=False)
            pd.DataFrame(updated_readouts).to_csv(old_exp, index=False)
    
    if VERBOSE: print("-SUCCESS- 2. Record files updated with CDD-Vault information.")
    
    # obtain list of ligands that need to be docked 
    protein_name = os.path.splittext(os.path.basename(pdb))[0]
    print("3. Checking exisiting predictions and performing Docking with Boltz-2 Covalent.")
    pred_rec = os.path.join(record_path, 'predictions.csv')
    if os.path.exists(pred_rec):
        pred_df = pd.read_csv(pred_rec)
        no_pred_prot = update_predictions.fetch_new(pred_df, updated_metadata, protein_name)

        error_csv = os.path.join(record_path, 'errored.csv')
        if os.path.exists(error_csv): # update list to disinclude errored ligands 
            error_df = pd.read_csv(error_csv)
            no_pred_prot = update_predictions.check_attempted(error_df, no_pred_prot) 
        else: 
            if VERBOSE: print("No errored compounds. Attempting to Dock all compounds.")

    elif not os.path.exists(pred_rec):
        # simply assume no prediction was ever made 
        if VERBOSE: print("predictions.csv was not found. Attempting to Dock all compounds.")
        dock_compounds = updated_metadata[['substance_id', 'inchi_key', 'smiles']]

        # call submit job  
        submit_job.run_boltz_cov(prot_file=pdb, ligand_df=dock_compounds, boltz_cache=boltz_cache, 
                      res_idx=res_idx, ligand_chain=ligand_chain, VERBOSE=VERBOSE, 
                      output_dir=output_dir, slurm_template=slurm_template, msa_path=msa_path)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_file", required=True)
    args = parser.parse_args()

    main(args)