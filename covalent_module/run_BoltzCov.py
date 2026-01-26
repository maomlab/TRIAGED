import os
import pandas as pd
import sys
import json 
import argparse
from BoltzCov.update_ligands import pull_cdd_ligs

def read_json_args(json_file):
    '''
    Reads the json input file with arguments required to dock covalently with Boltz-2.
    
    :param json_file: Path to json file.
    '''

    with open(json_file, 'r') as jf:
        run_arguments = json.load(jf)

        record_csv = os.path.expandvars(run_arguments["RECORD_CSV"]) # can be empty for first time run since ligands are pulled from cdd but path must be given 

        cdd_api_key = run_arguments.get("CDD_API_KEY", None)
        vault_id = run_arguments.get("VAULT_ID", None) # number of the vault in CDD not the compound ids in CDD
        readout_query = run_arguments.get("READOUT_KEY", None) # is a dict 
        mol_query = run_arguments.get("MOL_QUERY", None) # is a dict 
        readout_id = run_arguments.get("READOUT_ID", None) # is the experiment ID (either tgcpl 1125023 or hscpl 1125024 here)
        
        pdb = os.path.expandvars(run_arguments["PDB"])
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
    record_csv,
    (cdd_api_key, vault_id, readout_query, mol_query, readout_id),
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
    (record_csv,
    (cdd_api_key, vault_id, readout_query, mol_query, readout_id),
    (pdb, res_idx, ligand_chain, msa_path, boltz_cache, slurm_template),
    (VERBOSE, output_dir)
    ) = read_json_args(args.json_file)

    if any(x is None for x in (
    record_csv, pdb, res_idx, ligand_chain, boltz_cache, slurm_template,
    cdd_api_key, vault_id, readout_query, mol_query, readout_id, output_dir
    )):
        raise ValueError("Please make sure all required arguments are given in the input JSON.")

    # pull ligand information and experiment readouts 
    if VERBOSE: 
        print("Pulling ligands from CDD vault using the following queries:" \
            f"{readout_query}" \
            f"{mol_query}")
        
    readouts, molecules = pull_cdd_ligs.cdd_query(API_KEY=cdd_api_key, 
                                    VAULT_ID=7171, 
                                    readout_query=readout_query, 
                                    mol_query=mol_query)
    
    readouts_and_ids = pull_cdd_ligs.reorg_query(readouts, molecules, filter=False)

    compound_rec = pd.read_csv(args.record_csv)
    for idx, _ in compound_rec.iterrows():
        compound_rec.at[idx, 'comments'] = '' # clear all comments 

    # add cdd vault exp values 
    exp_update = pull_cdd_ligs.add_experiment_reads(cdd_query=readouts_and_ids, compound_rec_df=compound_rec)
    # add smiles from cdd vault
    smiles_update = pull_cdd_ligs.add_smiles(cdd_query=readouts_and_ids, compound_rec_df=exp_update)
    # add boltz predictions
    for pred_dir in args.pred_dirs:
        dir_name = os.path.basename(os.path.normpath(pred_dir))
        if 'tgcpl' in dir_name:
            smiles_update = pull_cdd_ligs.add_predictions(predictions_dir=pred_dir, compound_rec_df=smiles_update, target='tgcpl')
        elif 'hscpl' in dir_name:
            smiles_update = pull_cdd_ligs.add_predictions(predictions_dir=pred_dir, compound_rec_df=smiles_update, target='hscpl')
    
    smiles_update.to_csv(args.output)

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--json_file", required=True)
    args = parser.parse_args()

    main(args)